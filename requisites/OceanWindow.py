# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'OceanWindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1650, 809)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.morewidget = Ocean(self.centralwidget)
        self.morewidget.setGeometry(QtCore.QRect(10, 10, 1051, 741))
        self.morewidget.setObjectName("morewidget")
        self.Start_Button = QtWidgets.QPushButton(self.centralwidget)
        self.Start_Button.setGeometry(QtCore.QRect(1090, 680, 541, 71))
        font = QtGui.QFont()
        font.setPointSize(27)
        self.Start_Button.setFont(font)
        self.Start_Button.setObjectName("Start_Button")
        self.ForceSlider = QtWidgets.QSlider(self.centralwidget)
        self.ForceSlider.setGeometry(QtCore.QRect(1160, 490, 22, 160))
        self.ForceSlider.setMaximum(100)
        self.ForceSlider.setOrientation(QtCore.Qt.Vertical)
        self.ForceSlider.setObjectName("ForceSlider")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(1110, 460, 151, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(1110, 650, 131, 21))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.alpha_dial = QtWidgets.QDial(self.centralwidget)
        self.alpha_dial.setGeometry(QtCore.QRect(1380, 500, 151, 131))
        self.alpha_dial.setMinimum(-157)
        self.alpha_dial.setMaximum(157)
        self.alpha_dial.setObjectName("alpha_dial")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(1280, 550, 111, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(1530, 550, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(1070, 0, 211, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.m_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.m_Edit.setGeometry(QtCore.QRect(1100, 40, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.m_Edit.setFont(font)
        self.m_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.m_Edit.setObjectName("m_Edit")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(1070, 40, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.a_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.a_Edit.setGeometry(QtCore.QRect(1250, 40, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.a_Edit.setFont(font)
        self.a_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.a_Edit.setObjectName("a_Edit")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(1220, 40, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.h_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.h_Edit.setGeometry(QtCore.QRect(1400, 40, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.h_Edit.setFont(font)
        self.h_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.h_Edit.setObjectName("h_Edit")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(1370, 40, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.J_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.J_Edit.setGeometry(QtCore.QRect(1550, 40, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.J_Edit.setFont(font)
        self.J_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.J_Edit.setObjectName("J_Edit")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(1520, 40, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(1070, 80, 301, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.Cl_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Cl_Edit.setGeometry(QtCore.QRect(1100, 120, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Cl_Edit.setFont(font)
        self.Cl_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Cl_Edit.setObjectName("Cl_Edit")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(1070, 120, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.Cb_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Cb_Edit.setGeometry(QtCore.QRect(1250, 120, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Cb_Edit.setFont(font)
        self.Cb_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Cb_Edit.setObjectName("Cb_Edit")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(1220, 120, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(1360, 120, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.Cvr_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Cvr_Edit.setGeometry(QtCore.QRect(1400, 120, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Cvr_Edit.setFont(font)
        self.Cvr_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Cvr_Edit.setObjectName("Cvr_Edit")
        self.Rho_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Rho_Edit.setGeometry(QtCore.QRect(1550, 120, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Rho_Edit.setFont(font)
        self.Rho_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Rho_Edit.setObjectName("Rho_Edit")
        self.label_14 = QtWidgets.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(1510, 120, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(1520, 80, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(1370, 80, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.Sb_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Sb_Edit.setGeometry(QtCore.QRect(1550, 80, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Sb_Edit.setFont(font)
        self.Sb_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Sb_Edit.setObjectName("Sb_Edit")
        self.Sl_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Sl_Edit.setGeometry(QtCore.QRect(1400, 80, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Sl_Edit.setFont(font)
        self.Sl_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Sl_Edit.setObjectName("Sl_Edit")
        self.X_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.X_Edit.setGeometry(QtCore.QRect(1100, 370, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.X_Edit.setFont(font)
        self.X_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.X_Edit.setObjectName("X_Edit")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(1070, 370, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.Y_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Y_Edit.setGeometry(QtCore.QRect(1190, 370, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Y_Edit.setFont(font)
        self.Y_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Y_Edit.setObjectName("Y_Edit")
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(1160, 370, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_18.setFont(font)
        self.label_18.setObjectName("label_18")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(1240, 370, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_19.setFont(font)
        self.label_19.setObjectName("label_19")
        self.Phi_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Phi_Edit.setGeometry(QtCore.QRect(1280, 370, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Phi_Edit.setFont(font)
        self.Phi_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Phi_Edit.setObjectName("Phi_Edit")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        self.label_20.setGeometry(QtCore.QRect(1070, 330, 311, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.Vy_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Vy_Edit.setGeometry(QtCore.QRect(1190, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Vy_Edit.setFont(font)
        self.Vy_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Vy_Edit.setObjectName("Vy_Edit")
        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        self.label_21.setGeometry(QtCore.QRect(1250, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.Vx_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Vx_Edit.setGeometry(QtCore.QRect(1100, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Vx_Edit.setFont(font)
        self.Vx_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Vx_Edit.setObjectName("Vx_Edit")
        self.label_22 = QtWidgets.QLabel(self.centralwidget)
        self.label_22.setGeometry(QtCore.QRect(1160, 410, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_22.setFont(font)
        self.label_22.setObjectName("label_22")
        self.Omega_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.Omega_Edit.setGeometry(QtCore.QRect(1280, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Omega_Edit.setFont(font)
        self.Omega_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.Omega_Edit.setObjectName("Omega_Edit")
        self.label_23 = QtWidgets.QLabel(self.centralwidget)
        self.label_23.setGeometry(QtCore.QRect(1070, 410, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.label_28 = QtWidgets.QLabel(self.centralwidget)
        self.label_28.setGeometry(QtCore.QRect(1350, 330, 311, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_28.setFont(font)
        self.label_28.setObjectName("label_28")
        self.OmegaT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.OmegaT_Edit.setGeometry(QtCore.QRect(1580, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.OmegaT_Edit.setFont(font)
        self.OmegaT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.OmegaT_Edit.setObjectName("OmegaT_Edit")
        self.PhiT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.PhiT_Edit.setGeometry(QtCore.QRect(1580, 370, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.PhiT_Edit.setFont(font)
        self.PhiT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.PhiT_Edit.setObjectName("PhiT_Edit")
        self.label_24 = QtWidgets.QLabel(self.centralwidget)
        self.label_24.setGeometry(QtCore.QRect(1550, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_24.setFont(font)
        self.label_24.setObjectName("label_24")
        self.XT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.XT_Edit.setGeometry(QtCore.QRect(1400, 370, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.XT_Edit.setFont(font)
        self.XT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.XT_Edit.setObjectName("XT_Edit")
        self.label_25 = QtWidgets.QLabel(self.centralwidget)
        self.label_25.setGeometry(QtCore.QRect(1460, 370, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_25.setFont(font)
        self.label_25.setObjectName("label_25")
        self.label_26 = QtWidgets.QLabel(self.centralwidget)
        self.label_26.setGeometry(QtCore.QRect(1540, 370, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_26.setFont(font)
        self.label_26.setObjectName("label_26")
        self.VxT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.VxT_Edit.setGeometry(QtCore.QRect(1400, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.VxT_Edit.setFont(font)
        self.VxT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.VxT_Edit.setObjectName("VxT_Edit")
        self.label_27 = QtWidgets.QLabel(self.centralwidget)
        self.label_27.setGeometry(QtCore.QRect(1370, 370, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_27.setFont(font)
        self.label_27.setObjectName("label_27")
        self.YT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.YT_Edit.setGeometry(QtCore.QRect(1490, 370, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.YT_Edit.setFont(font)
        self.YT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.YT_Edit.setObjectName("YT_Edit")
        self.label_29 = QtWidgets.QLabel(self.centralwidget)
        self.label_29.setGeometry(QtCore.QRect(1460, 410, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_29.setFont(font)
        self.label_29.setObjectName("label_29")
        self.VyT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.VyT_Edit.setGeometry(QtCore.QRect(1490, 410, 51, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.VyT_Edit.setFont(font)
        self.VyT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.VyT_Edit.setObjectName("VyT_Edit")
        self.label_30 = QtWidgets.QLabel(self.centralwidget)
        self.label_30.setGeometry(QtCore.QRect(1370, 410, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_30.setFont(font)
        self.label_30.setObjectName("label_30")
        self.SlT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.SlT_Edit.setGeometry(QtCore.QRect(1400, 240, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.SlT_Edit.setFont(font)
        self.SlT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.SlT_Edit.setObjectName("SlT_Edit")
        self.mT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.mT_Edit.setGeometry(QtCore.QRect(1100, 200, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.mT_Edit.setFont(font)
        self.mT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.mT_Edit.setObjectName("mT_Edit")
        self.label_31 = QtWidgets.QLabel(self.centralwidget)
        self.label_31.setGeometry(QtCore.QRect(1220, 200, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_31.setFont(font)
        self.label_31.setObjectName("label_31")
        self.ClT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.ClT_Edit.setGeometry(QtCore.QRect(1100, 280, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.ClT_Edit.setFont(font)
        self.ClT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.ClT_Edit.setObjectName("ClT_Edit")
        self.label_32 = QtWidgets.QLabel(self.centralwidget)
        self.label_32.setGeometry(QtCore.QRect(1070, 160, 251, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_32.setFont(font)
        self.label_32.setObjectName("label_32")
        self.label_33 = QtWidgets.QLabel(self.centralwidget)
        self.label_33.setGeometry(QtCore.QRect(1070, 280, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_33.setFont(font)
        self.label_33.setObjectName("label_33")
        self.VoT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.VoT_Edit.setGeometry(QtCore.QRect(1550, 280, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.VoT_Edit.setFont(font)
        self.VoT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.VoT_Edit.setObjectName("VoT_Edit")
        self.label_34 = QtWidgets.QLabel(self.centralwidget)
        self.label_34.setGeometry(QtCore.QRect(1520, 200, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_34.setFont(font)
        self.label_34.setObjectName("label_34")
        self.RT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.RT_Edit.setGeometry(QtCore.QRect(1400, 200, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.RT_Edit.setFont(font)
        self.RT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.RT_Edit.setObjectName("RT_Edit")
        self.label_35 = QtWidgets.QLabel(self.centralwidget)
        self.label_35.setGeometry(QtCore.QRect(1220, 280, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_35.setFont(font)
        self.label_35.setObjectName("label_35")
        self.label_36 = QtWidgets.QLabel(self.centralwidget)
        self.label_36.setGeometry(QtCore.QRect(1360, 280, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_36.setFont(font)
        self.label_36.setObjectName("label_36")
        self.CvrT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.CvrT_Edit.setGeometry(QtCore.QRect(1400, 280, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.CvrT_Edit.setFont(font)
        self.CvrT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.CvrT_Edit.setObjectName("CvrT_Edit")
        self.label_37 = QtWidgets.QLabel(self.centralwidget)
        self.label_37.setGeometry(QtCore.QRect(1070, 240, 301, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_37.setFont(font)
        self.label_37.setObjectName("label_37")
        self.label_38 = QtWidgets.QLabel(self.centralwidget)
        self.label_38.setGeometry(QtCore.QRect(1370, 240, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_38.setFont(font)
        self.label_38.setObjectName("label_38")
        self.label_39 = QtWidgets.QLabel(self.centralwidget)
        self.label_39.setGeometry(QtCore.QRect(1370, 200, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_39.setFont(font)
        self.label_39.setObjectName("label_39")
        self.label_40 = QtWidgets.QLabel(self.centralwidget)
        self.label_40.setGeometry(QtCore.QRect(1520, 240, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_40.setFont(font)
        self.label_40.setObjectName("label_40")
        self.label_41 = QtWidgets.QLabel(self.centralwidget)
        self.label_41.setGeometry(QtCore.QRect(1510, 280, 41, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_41.setFont(font)
        self.label_41.setObjectName("label_41")
        self.aT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.aT_Edit.setGeometry(QtCore.QRect(1250, 200, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.aT_Edit.setFont(font)
        self.aT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.aT_Edit.setObjectName("aT_Edit")
        self.SbT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.SbT_Edit.setGeometry(QtCore.QRect(1550, 240, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.SbT_Edit.setFont(font)
        self.SbT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.SbT_Edit.setObjectName("SbT_Edit")
        self.CbT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.CbT_Edit.setGeometry(QtCore.QRect(1250, 280, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.CbT_Edit.setFont(font)
        self.CbT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.CbT_Edit.setObjectName("CbT_Edit")
        self.label_42 = QtWidgets.QLabel(self.centralwidget)
        self.label_42.setGeometry(QtCore.QRect(1070, 200, 31, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_42.setFont(font)
        self.label_42.setObjectName("label_42")
        self.JT_Edit = QtWidgets.QLineEdit(self.centralwidget)
        self.JT_Edit.setGeometry(QtCore.QRect(1550, 200, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.JT_Edit.setFont(font)
        self.JT_Edit.setAlignment(QtCore.Qt.AlignCenter)
        self.JT_Edit.setObjectName("JT_Edit")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1650, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.Start_Button.setText(_translate("MainWindow", "Поехали!"))
        self.label.setText(_translate("MainWindow", "Полный вперёд"))
        self.label_2.setText(_translate("MainWindow", "Стоп машина"))
        self.label_3.setText(_translate("MainWindow", "Право руля"))
        self.label_4.setText(_translate("MainWindow", "Лево руля"))
        self.label_5.setText(_translate("MainWindow", "Параметры лодки:"))
        self.m_Edit.setText(_translate("MainWindow", "400"))
        self.label_6.setText(_translate("MainWindow", "m:"))
        self.a_Edit.setText(_translate("MainWindow", "1.5"))
        self.label_7.setText(_translate("MainWindow", "a:"))
        self.h_Edit.setText(_translate("MainWindow", "0.3"))
        self.label_8.setText(_translate("MainWindow", "h:"))
        self.J_Edit.setText(_translate("MainWindow", "530"))
        self.label_9.setText(_translate("MainWindow", "J:"))
        self.label_10.setText(_translate("MainWindow", "Параметры сопротивления:"))
        self.Cl_Edit.setText(_translate("MainWindow", "0.2"))
        self.label_11.setText(_translate("MainWindow", "Cl:"))
        self.Cb_Edit.setText(_translate("MainWindow", "0.9"))
        self.label_12.setText(_translate("MainWindow", "Cb:"))
        self.label_13.setText(_translate("MainWindow", "Cvr:"))
        self.Cvr_Edit.setText(_translate("MainWindow", "2"))
        self.Rho_Edit.setText(_translate("MainWindow", "1000"))
        self.label_14.setText(_translate("MainWindow", "Rho:"))
        self.label_15.setText(_translate("MainWindow", "Sb:"))
        self.label_16.setText(_translate("MainWindow", "Sl:"))
        self.Sb_Edit.setText(_translate("MainWindow", "1.8"))
        self.Sl_Edit.setText(_translate("MainWindow", "0.15"))
        self.X_Edit.setText(_translate("MainWindow", "0"))
        self.label_17.setText(_translate("MainWindow", "X:"))
        self.Y_Edit.setText(_translate("MainWindow", "0"))
        self.label_18.setText(_translate("MainWindow", "Y:"))
        self.label_19.setText(_translate("MainWindow", "Phi:"))
        self.Phi_Edit.setText(_translate("MainWindow", "0"))
        self.label_20.setText(_translate("MainWindow", "Начальное состояние лодки:"))
        self.Vy_Edit.setText(_translate("MainWindow", "0"))
        self.label_21.setText(_translate("MainWindow", " ω:"))
        self.Vx_Edit.setText(_translate("MainWindow", "0"))
        self.label_22.setText(_translate("MainWindow", "Vy:"))
        self.Omega_Edit.setText(_translate("MainWindow", "0"))
        self.label_23.setText(_translate("MainWindow", "Vx:"))
        self.label_28.setText(_translate("MainWindow", "Начальное состояние торпеды:"))
        self.OmegaT_Edit.setText(_translate("MainWindow", "0"))
        self.PhiT_Edit.setText(_translate("MainWindow", "0"))
        self.label_24.setText(_translate("MainWindow", " ω:"))
        self.XT_Edit.setText(_translate("MainWindow", "-60"))
        self.label_25.setText(_translate("MainWindow", "Y:"))
        self.label_26.setText(_translate("MainWindow", "Phi:"))
        self.VxT_Edit.setText(_translate("MainWindow", "0"))
        self.label_27.setText(_translate("MainWindow", "X:"))
        self.YT_Edit.setText(_translate("MainWindow", "0"))
        self.label_29.setText(_translate("MainWindow", "Vy:"))
        self.VyT_Edit.setText(_translate("MainWindow", "0"))
        self.label_30.setText(_translate("MainWindow", "Vx:"))
        self.SlT_Edit.setText(_translate("MainWindow", "0.15"))
        self.mT_Edit.setText(_translate("MainWindow", "400"))
        self.label_31.setText(_translate("MainWindow", "a:"))
        self.ClT_Edit.setText(_translate("MainWindow", "0.2"))
        self.label_32.setText(_translate("MainWindow", "Параметры торпеды:"))
        self.label_33.setText(_translate("MainWindow", "Cl:"))
        self.VoT_Edit.setText(_translate("MainWindow", "4"))
        self.label_34.setText(_translate("MainWindow", "J:"))
        self.RT_Edit.setText(_translate("MainWindow", "0.3"))
        self.label_35.setText(_translate("MainWindow", "Cb:"))
        self.label_36.setText(_translate("MainWindow", "Cvr:"))
        self.CvrT_Edit.setText(_translate("MainWindow", "2"))
        self.label_37.setText(_translate("MainWindow", "Параметры сопротивления:"))
        self.label_38.setText(_translate("MainWindow", "Sl:"))
        self.label_39.setText(_translate("MainWindow", "R:"))
        self.label_40.setText(_translate("MainWindow", "Sb:"))
        self.label_41.setText(_translate("MainWindow", "VoT:"))
        self.aT_Edit.setText(_translate("MainWindow", "1.5"))
        self.SbT_Edit.setText(_translate("MainWindow", "1.8"))
        self.CbT_Edit.setText(_translate("MainWindow", "0.9"))
        self.label_42.setText(_translate("MainWindow", "m:"))
        self.JT_Edit.setText(_translate("MainWindow", "530"))
from ocean import Ocean
