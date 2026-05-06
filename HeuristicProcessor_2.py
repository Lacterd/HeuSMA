# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 20:24:43 2023

@author: lacter
"""
import pandas as pd
import time
from progress.bar import Bar
import pyopenms
import math
import numpy as np
import bisect
import re
import os
from scipy.optimize import linear_sum_assignment
from scipy.spatial import KDTree
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder, StandardScaler
from scipy.stats import chi2
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
#from sklearn.linear_model import ElasticNetCV
from sklearn.cross_decomposition import PLSRegression
#from scipy.cluster.hierarchy import linkage, leaves_list
#import plotly.figure_factory as ff
#from plotly.subplots import make_subplots
#from sklearn.model_selection import train_test_split
#from sklearn.metrics import accuracy_score
#from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import warnings
import glob
#import json
import sqlite3
'''---------'''
import sys
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QGridLayout, QLabel, QTableWidget, QHeaderView,
    QPushButton, QLineEdit, QProgressBar, QDesktopWidget, QFileDialog,
    QTableWidgetItem, QApplication, QComboBox, QDialog
)
from PyQt5.QtGui import QFont,QIcon
from PyQt5.QtCore import QThread, pyqtSignal
from concurrent.futures import ProcessPoolExecutor, as_completed
#from PyQt5.QtCore import Qt
from PyQt5 import QtWidgets
#from PyQt5 import QtCore
#import pathos
import multiprocessing as mp
import pickle
import plotly.graph_objects as go
import plotly.express as px
#import plotly.io as pio
import networkx as nx
#from collections import defaultdict

class HeuristicProcessorUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
        # 设置标题
        self.resize(700,900)
        self.centerWindow()
        
        # 设置参数字典
        self.Data_params={}
        self.Calculate_params={'MS1_Tor':0.000010,'smooth':5,'min_Int':10000,'Points':17}
    
    def CloseEvent(self,event):
        os._exit(0)

    def initUI(self):
        globallayout = QVBoxLayout()
        # Data import
        Data_import_Widget = QWidget()
        Data_import_Layout = QGridLayout()
        # Union
        self.Label_SampleSelect = QLabel('Select data')
        self.Label_RIISSelect = QLabel('Select Calibrants data')
        self.Label_HeuristicListSelect = QLabel('Select heuristic list')
        self.Label_Params_Select = QLabel('Set params')
        self.Label_MS2_Params = QLabel('Annotation and Cluster')
        self.TextBrowser_SampleSelect = QTableWidget()
        self.TextBrowser_SampleSelect.setColumnCount(4)
        self.TextBrowser_SampleSelect.setHorizontalHeaderLabels(['Name','Type','Group','Index'])
        self.TextBrowser_SampleSelect.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.TextBrowser_SampleSelect.horizontalHeader().setStretchLastSection(True)
        self.TextBrowser_RIISSelect = QLabel('*')
        self.TextBrowser_HeuristicListSelect = QLabel('*')
        self.PushButton_SampleSelect = QPushButton('Select')
        self.PushButton_RIISSelect = QPushButton('Select')
        self.PushButton_HeuristicListSelect = QPushButton('Select')
        self.PushButton_SampleSelect.setToolTip('Select sample *.mzML documents')
        self.PushButton_RIISSelect.setToolTip('Select RIIS *.xlsx or *.xls document')
        self.PushButton_HeuristicListSelect.setToolTip('Select Heuristic List *.xlsx or *.xls document')
        # Slot
        self.PushButton_SampleSelect.clicked.connect(self.SelectSampleFile)
        self.PushButton_RIISSelect.clicked.connect(self.SelectRIISFile)
        self.PushButton_HeuristicListSelect.clicked.connect(self.SelectHeuristicList)
        Data_import_Layout.addWidget(self.Label_SampleSelect,1,0)
        Data_import_Layout.addWidget(self.TextBrowser_SampleSelect, 2, 0) # row 1， column 0
        Data_import_Layout.addWidget(self.PushButton_SampleSelect, 2, 1) # row 1， column 1
        Data_import_Layout.addWidget(self.Label_RIISSelect,5,0)
        Data_import_Layout.addWidget(self.TextBrowser_RIISSelect, 6, 0)
        Data_import_Layout.addWidget(self.PushButton_RIISSelect, 6, 1) # 行，列，行高，列宽
        Data_import_Layout.addWidget(self.Label_HeuristicListSelect,7,0)
        Data_import_Layout.addWidget(self.TextBrowser_HeuristicListSelect, 8, 0)
        Data_import_Layout.addWidget(self.PushButton_HeuristicListSelect, 8, 1) # 行，列，行高，列宽
        Data_import_Widget.setLayout(Data_import_Layout)
        # Params setting
        Params_setting_Widget = QWidget()
        Params_setting_Layout = QGridLayout()
        # Union
        self.Lable_MS_Tor = QLabel('MS1 Tolerance(ppm)')
        self.Lable_RI_Tor = QLabel('RI Tolerance(%)')
        self.Lable_Int_min = QLabel('Min intensity')
        self.Lable_Point = QLabel('Min scan points')
        self.Lable_SN = QLabel('Signal/Noise')
        self.Lable_SB = QLabel('Signal/Blank')
        self.Lable_RT_Tor = QLabel('RT Tolerance(min)')
        self.LineEdit_MS_Tor = QLineEdit('10')
        self.LineEdit_RI_Tor = QLineEdit('2')
        self.LineEdit_Int_min = QLineEdit('10000')
        self.LineEdit_Point = QLineEdit('17')
        self.LineEdit_SN = QLineEdit('10')
        self.LineEdit_SB = QLineEdit('3')
        self.LineEdit_RT_Tor = QLineEdit('0.1')
        self.LineEdit_MS_Tor.setToolTip('Under this threshold will be recognized as same m/z')
        self.LineEdit_RI_Tor.setToolTip('Under this threshold will be recognized as same RI')
        self.LineEdit_Int_min.setToolTip('Under this threshold will be deleted')
        self.LineEdit_Point.setToolTip('Scan points under this threshold will be deleted')
        self.LineEdit_SN.setToolTip('Under this threshold will be deleted')
        self.LineEdit_SB.setToolTip('Under this threshold will be deleted')
        self.LineEdit_RT_Tor.setToolTip('Under this threshold will be recognized as same RT')
        Params_setting_Layout.addWidget(self.Label_Params_Select,0,0)
        #Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.LineEdit_MS_Tor,1,1)
        Params_setting_Layout.addWidget(self.Lable_RI_Tor,1,2)
        Params_setting_Layout.addWidget(self.LineEdit_RI_Tor,1,3)
        Params_setting_Layout.addWidget(self.Lable_Int_min,1,4)
        Params_setting_Layout.addWidget(self.LineEdit_Int_min,1,5)
        Params_setting_Layout.addWidget(self.Lable_Point,2,0)
        Params_setting_Layout.addWidget(self.LineEdit_Point,2,1)
        Params_setting_Layout.addWidget(self.Lable_SN,2,2)
        Params_setting_Layout.addWidget(self.LineEdit_SN,2,3)
        Params_setting_Layout.addWidget(self.Lable_SB,2,4)
        Params_setting_Layout.addWidget(self.LineEdit_SB,2,5)
        Params_setting_Layout.addWidget(self.Lable_RT_Tor,3,0)
        Params_setting_Layout.addWidget(self.LineEdit_RT_Tor,3,1)
        Params_setting_Widget.setLayout(Params_setting_Layout)
        # Annotation
        MS2_setting_Widget = QWidget()
        MS2_setting_Layout = QGridLayout()
        # Union
        self.Lable_Annotation_MS2_Number = QLabel('Max match ions')
        self.Lable_min_cos = QLabel('Min pair Cos')
        #self.Lable_max_K = QLabel('Max network points')
        self.Lable_DataBase = QLabel('Database')
        self.LineEdit_Annotation_MS2_Number = QLineEdit('5')
        self.LineEdit_Annotation_MS2_Number.setToolTip('At most this number of ions will be used to calculate')
        self.LineEdit_min_cos = QLineEdit('0.7')
        self.LineEdit_min_cos.setToolTip('Under this threshold will not be matched')
        #self.LineEdit_max_K = QLineEdit('20')
        self.PushButton_DB_Select = QPushButton('Select')
        self.PushButton_DB_Select.setToolTip('Select database which will be used')
        self.TextBrowser_DB_Select = QLabel('*')
        self.PushButton_DB_Select.clicked.connect(self.SelectDBFile)
        MS2_setting_Layout.addWidget(self.Label_MS2_Params,0,0)
        MS2_setting_Layout.addWidget(self.Lable_Annotation_MS2_Number,1,0)
        MS2_setting_Layout.addWidget(self.LineEdit_Annotation_MS2_Number,1,1)
        MS2_setting_Layout.addWidget(self.Lable_min_cos,1,2)
        MS2_setting_Layout.addWidget(self.LineEdit_min_cos,1,3)
        #MS2_setting_Layout.addWidget(self.Lable_max_K,1,4)
        #MS2_setting_Layout.addWidget(self.LineEdit_max_K,1,5)
        MS2_setting_Layout.addWidget(self.Lable_DataBase,2,0)
        MS2_setting_Layout.addWidget(self.TextBrowser_DB_Select,2,1)
        MS2_setting_Layout.addWidget(self.PushButton_DB_Select,2,3)
        MS2_setting_Widget.setLayout(MS2_setting_Layout)
        # Run button
        Run_button_Widget = QWidget()
        Run_button_Layout = QGridLayout()
        # Union
        self.PushButton_Run = QPushButton('Run')
        self.Label_Run1 = QLabel('')
        self.Label_Run2 = QLabel('')
        self.Label_Run3 = QLabel('')
        Run_button_Layout.addWidget(self.Label_Run1,0,0)
        Run_button_Layout.addWidget(self.Label_Run2,0,1)
        Run_button_Layout.addWidget(self.Label_Run3,0,2)
        Run_button_Layout.addWidget(self.PushButton_Run,0,3)
        Run_button_Widget.setLayout(Run_button_Layout)
        self.PushButton_Run.clicked.connect(self.Run)
        
        # statusbar 
        StatusBar_Widget = QWidget()
        StatusBar_Layout = QGridLayout()
        # Union
        
        self.Label_process = QLabel('Processing bar')
        self.Label_process_sub = QLabel('Step')
        self.process_bar = QProgressBar()
        self.process_bar.setStyleSheet("QProgressBar { border: 2px solid grey; border-radius: 5px; color: rgb(20,20,20);  background-color: #FFFFFF; text-align: center;}QProgressBar::chunk {background-color: rgb(100,200,200); border-radius: 10px; margin: 0.1px;  width: 1px;}")
        font = QFont()
        font.setBold(True)
        font.setWeight(30)
        self.process_bar.setFont(font)
        self.process_bar.setMaximum(100)
        self.process_bar.setMinimum(0)
        self.process_bar.setValue(0)
        StatusBar_Layout.addWidget(self.Label_process,0,0)
        StatusBar_Layout.addWidget(self.Label_process_sub,1,0)
        StatusBar_Layout.addWidget(self.process_bar,1,1)
        StatusBar_Widget.setLayout(StatusBar_Layout)
        # Global Layout setting
        globallayout.addWidget(Data_import_Widget)
        globallayout.addWidget(Params_setting_Widget)
        globallayout.addWidget(MS2_setting_Widget)
        globallayout.addWidget(Run_button_Widget)
        globallayout.addWidget(StatusBar_Widget)
        self.setLayout(globallayout)
        self.setWindowTitle('Heuristic Processor')
        
    def centerWindow(self):
        screen = QDesktopWidget().screenGeometry()
        size  = self.geometry()
        LeftValue  = int((screen.width()-size.width())/2)
        TopValue = int((screen.height()-size.height())/2)
        self.move(LeftValue,TopValue)
    
    def SelectSampleFile(self):
        FileName,FileType = QFileDialog.getOpenFileNames(self,"Select files",os.getcwd(),"mzML Files (*.mzML);;Excel File (*.xlsx)")
        if len(FileName) > 0 and FileName[0].endswith('.mzML'):
            for i in FileName:
                if i.endswith('.mzML'):
                    Combox_Type = QComboBox()
                    Combox_Type.addItems(['Sample','QC','Blank','MS2'])
                    Combox_Type.setCurrentText('Sample')
                    Combox_Type.currentIndexChanged.connect(self.SampleTypeChange)
                    LineEdit_Index= QLineEdit(str(self.TextBrowser_SampleSelect.rowCount()+1))
                    LineEdit_Index.textChanged.connect(self.SampleIndexChange)
                    LineEdit_Group= QLineEdit('1')
                    LineEdit_Group.textChanged.connect(self.SampleGroupChange)
                    name_begin = i.rfind('/')
                    if name_begin == -1:
                        name_begin = i.rfind('\\')
                    name_end = len(i)
                    self.Data_params[i[name_begin+1:name_end-5]] = {'Name':i[name_begin+1:name_end-5],'Path':i,'Type':'Sample','Group':1,'Index':str(self.TextBrowser_SampleSelect.rowCount()),'Combox':Combox_Type,'LineEdit_Index':LineEdit_Index,'LineEdit_Group':LineEdit_Group,'rowCount':self.TextBrowser_SampleSelect.rowCount()}
                    self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
                    self.TextBrowser_SampleSelect.setItem(self.TextBrowser_SampleSelect.rowCount()-1,0,QTableWidgetItem(i[name_begin+1:name_end-5]))
                    self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,self.Data_params[i[name_begin+1:name_end-5]]['Combox'])
                    self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,2,self.Data_params[i[name_begin+1:name_end-5]]['LineEdit_Group'])
                    self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,3,self.Data_params[i[name_begin+1:name_end-5]]['LineEdit_Index'])
                    QApplication.processEvents()
        elif len(FileName) == 1 and FileName[0].lower().endswith('.xlsx'):
            SampleData_excel = pd.read_excel(FileName[0])
            for i in range(len(SampleData_excel)):
                Combox_Type = QComboBox()
                Combox_Type.addItems(['Sample','QC','Blank','MS2'])
                if SampleData_excel.loc[i,'SampleType'] == 'Sample':
                    Combox_Type.setCurrentText('Sample')
                elif SampleData_excel.loc[i,'SampleType'] == 'QC':
                    Combox_Type.setCurrentText('QC')
                elif SampleData_excel.loc[i,'SampleType'] == 'Blank':
                    Combox_Type.setCurrentText('Blank')
                elif SampleData_excel.loc[i,'SampleType'] == 'MS2':
                    Combox_Type.setCurrentText('MS2')
                Combox_Type.currentIndexChanged.connect(self.SampleTypeChange)
                if SampleData_excel.loc[i,'SampleType'] == 'Sample' or SampleData_excel.loc[i,'SampleType'] == 'QC':
                    LineEdit_Index= QLineEdit(str(int(SampleData_excel.loc[i,'Index'])))
                    LineEdit_Index.textChanged.connect(self.SampleIndexChange)
                    if SampleData_excel.loc[i,'SampleType'] == 'Sample':
                        LineEdit_Group= QLineEdit(SampleData_excel.loc[i,'Group'])
                    else:
                        LineEdit_Group= QLineEdit('QC')
                    LineEdit_Group.textChanged.connect(self.SampleGroupChange)
                name_begin = SampleData_excel.loc[i,'Path'].rfind('/')
                if name_begin == -1:
                    name_begin = SampleData_excel.loc[i,'Path'].rfind('\\')
                name_end = len(SampleData_excel.loc[i,'Path'])
                if SampleData_excel.loc[i,'SampleType'] == 'Sample':
                    self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]] = {'Name':SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5],'Path':SampleData_excel.loc[i,'Path'],'Type':SampleData_excel.loc[i,'SampleType'],'Group':SampleData_excel.loc[i,'Group'],'Index':str(int(SampleData_excel.loc[i,'Index'])),'Combox':Combox_Type,'LineEdit_Index':LineEdit_Index,'LineEdit_Group':LineEdit_Group,'rowCount':self.TextBrowser_SampleSelect.rowCount()}
                elif SampleData_excel.loc[i,'SampleType'] == 'QC':
                    self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]] = {'Name':SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5],'Path':SampleData_excel.loc[i,'Path'],'Type':'QC','Group':'QC','Index':str(int(SampleData_excel.loc[i,'Index'])),'Combox':Combox_Type,'LineEdit_Index':LineEdit_Index,'LineEdit_Group':LineEdit_Group,'rowCount':self.TextBrowser_SampleSelect.rowCount()}
                elif SampleData_excel.loc[i,'SampleType'] == 'Blank':
                    LineEdit_Index= QLineEdit('0')
                    LineEdit_Index.textChanged.connect(self.SampleIndexChange)
                    LineEdit_Group= QLineEdit('Blank')
                    LineEdit_Group.textChanged.connect(self.SampleGroupChange)
                    self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]] = {'Name':SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5],'Path':SampleData_excel.loc[i,'Path'],'Type':'Blank','Group':'Blank','Index':str(0),'Combox':Combox_Type,'LineEdit_Index':LineEdit_Index,'LineEdit_Group':LineEdit_Group,'rowCount':self.TextBrowser_SampleSelect.rowCount()}
                elif SampleData_excel.loc[i,'SampleType'] == 'MS2':
                    LineEdit_Index= QLineEdit('0')
                    LineEdit_Index.textChanged.connect(self.SampleIndexChange)
                    LineEdit_Group= QLineEdit('MS2')
                    LineEdit_Group.textChanged.connect(self.SampleGroupChange)
                    self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]] = {'Name':SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5],'Path':SampleData_excel.loc[i,'Path'],'Type':'MS2','Group':'MS2','Index':str(0),'Combox':Combox_Type,'LineEdit_Index':LineEdit_Index,'LineEdit_Group':LineEdit_Group,'rowCount':self.TextBrowser_SampleSelect.rowCount()}
                self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
                self.TextBrowser_SampleSelect.setItem(self.TextBrowser_SampleSelect.rowCount()-1,0,QTableWidgetItem(SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]))
                self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]]['Combox'])
                self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,2,self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]]['LineEdit_Group'])
                self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,3,self.Data_params[SampleData_excel.loc[i,'Path'][name_begin+1:name_end-5]]['LineEdit_Index'])
                QApplication.processEvents()
                
    def SelectRIISFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"选取文件",os.getcwd(),"Excel Files(*.xlsx)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_RIISSelect.setText(FileName[name_begin+1:name_end])
        self.RIISPath = FileName
        
    def SelectDBFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"选取文件",os.getcwd(),"Excel Files(*.db)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_DB_Select.setText(FileName[name_begin+1:name_end])
        self.DB_Path = FileName
        
    def SelectHeuristicList(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"选取文件",os.getcwd(),"Pickle Files(*.pkl)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_HeuristicListSelect.setText(FileName[name_begin+1:name_end])
        self.HeuristicListPath = FileName
    
    def SampleTypeChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['Combox']:
                self.Data_params[i]['Type'] = self.Data_params[i]['Combox'].currentText()
                if self.Data_params[i]['Combox'].currentText() == 'QC':
                    self.Data_params[i]['Group'] = 'QC'
                    self.Data_params[i]['LineEdit_Group'] = QLineEdit('QC')
                    self.Data_params[i]['LineEdit_Group'].setReadOnly(True)
                    self.TextBrowser_SampleSelect.setCellWidget(self.Data_params[i]['rowCount'],2,self.Data_params[i]['LineEdit_Group'])
                elif self.Data_params[i]['Combox'].currentText() == 'MS2':
                    self.Data_params[i]['Group'] = 'MS2'
                    self.Data_params[i]['LineEdit_Group'] = QLineEdit('MS2')
                    self.Data_params[i]['LineEdit_Group'].setReadOnly(True)
                    self.TextBrowser_SampleSelect.setCellWidget(self.Data_params[i]['rowCount'],2,self.Data_params[i]['LineEdit_Group'])
                elif self.Data_params[i]['Combox'].currentText() == 'Blank':
                    self.Data_params[i]['Group'] = 'Blank'
                    self.Data_params[i]['LineEdit_Group'] = QLineEdit('Blank')
                    self.Data_params[i]['LineEdit_Group'].setReadOnly(True)
                    self.TextBrowser_SampleSelect.setCellWidget(self.Data_params[i]['rowCount'],2,self.Data_params[i]['LineEdit_Group'])
                else:
                    self.Data_params[i]['LineEdit_Group'] = QLineEdit('1')
                    self.TextBrowser_SampleSelect.setCellWidget(self.Data_params[i]['rowCount'],2,self.Data_params[i]['LineEdit_Group'])
            QApplication.processEvents()
        
    def SampleIndexChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['LineEdit_Index']:
               self.Data_params[i]['Index'] = self.Data_params[i]['LineEdit_Index'].text()
               
    def SampleGroupChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['LineEdit_Group']:
               self.Data_params[i]['Group'] = self.Data_params[i]['LineEdit_Group'].text()
    
    def processbar_fresh(self,*arg):
        self.process_bar.setValue(self.process_bar.value()+1)
        QApplication.processEvents()
        
    def update_progress(self, value):
        self.process_bar.setValue(value)
        QApplication.processEvents()
        
    def on_PeakDetectWorker_done(self,result_dict):
        self.Sample_Pool = result_dict
        for i in self.SampleData.keys():
            self.SampleData[i]['Final_Peak_Detect'] = result_dict[i][0].copy()
            self.SampleData[i]['RIIS'] = result_dict[i][1].copy()
        temp_key = list(self.SampleData.keys())[0]
        temp_EMZDP = EazyMZDataProcess(self.SampleData[temp_key]['Path'])
        self.Polarity = temp_EMZDP.get_param('Polarity')
        ''' --- Process QC Data ---'''
        self.Label_process_sub.setText('Load QC Data')
        self.process_bar.setMaximum(100)
        self.process_bar.setValue(0)
        QApplication.processEvents()
        self.worker_QC = PeakDetectWorker(self.QCData, self.params)
        self.worker_QC.progress.connect(self.update_progress)
        self.worker_QC.finished.connect(self.on_QCDetectWorker_done)
        self.worker_QC.start()
    
    def on_QCDetectWorker_done(self,result_dict):
        self.QC_Pool = result_dict
        for i in self.QCData.keys():
            self.QCData[i]['Final_Peak_Detect'] = result_dict[i][0].copy()
            self.QCData[i]['RIIS'] = result_dict[i][1].copy()
        ''' --- Load Blank Data ---'''
        self.Label_process_sub.setText('Load Blank Data')
        self.process_bar.setMaximum(len(self.BlankData.keys()))
        self.process_bar.setValue(0)
        QApplication.processEvents()
        for i in self.BlankData.keys():
            self.process_bar.setValue(self.process_bar.value()+1)
            QApplication.processEvents()
            self.BlankData[i]['Data'] = EazyMZDataProcess(self.BlankData[i]['Path'])
            self.BlankData[i]['Data'].set_RIIS(self.RIISPath)
            self.BlankData[i]['Data'].set_param('MS1_Tor',float(self.LineEdit_MS_Tor.text())/1000000)
            self.BlankData[i]['Data'].set_param('RI_Tor',float(self.LineEdit_RI_Tor.text())/100)
            self.BlankData[i]['Data'].set_param('RT_Tor',float(self.LineEdit_RT_Tor.text())*60)
            self.BlankData[i]['Data'].set_param('min_Int',int(self.LineEdit_Int_min.text()))
            self.BlankData[i]['RIIS'] = self.BlankData[i]['Data'].RIIS.copy()
        ''' --- Load MS2 Data ---'''
        #self.Label_process.setText('Load MS2 Data')
        self.Label_process_sub.setText('Load MS2 Data')
        self.process_bar.setMaximum(len(self.MS2_Data.keys()))
        self.process_bar.setValue(0)
        self.worker_MS2 = MS2Worker(self.MS2_Data, self.params)
        self.worker_MS2.progress.connect(self.update_progress)
        self.worker_MS2.finished.connect(self.on_load_MS2_done)
        self.worker_MS2.start()
    
    def on_load_MS2_done(self,result_dict):
        self.MS2_Pool = result_dict
        for i in self.MS2_Data.keys():
            self.MS2_Data[i]['MS2_Data'] = result_dict[i].copy()
        ''' --- Alignment ---'''
        #self.Label_process.setText('Data Alignment')
        self.Label_process_sub.setText('Aligning')
        QApplication.processEvents()
        self.Align = DataAlignment(self.Label_process_sub,self.process_bar)
        for i in self.SampleData.keys():
            self.Align.add_Data(self.SampleData[i],Tag='Sample')
        for i in self.QCData.keys():
            self.Align.add_Data(self.QCData[i],Tag='QC')
        for i in self.BlankData.keys():
            self.Align.add_Data(self.BlankData[i],Tag='Blank')
        self.Align.set_param('RI_Alignment',True)
        self.Align.set_param('MZ_Tor',float(self.LineEdit_MS_Tor.text())/1000000)
        self.Align.set_param('RI_Tor',float(self.LineEdit_RI_Tor.text())/100)
        self.Align.RenewRefList()
        ''' Filter '''
        for i in self.ClassList:
            self.Align.RefList['Fill Group '+str(i)] = 0
            for ii in range(len(self.Align.RefList)):
                Number_Class = 0
                for iii in self.SampleClass.keys():
                    if self.SampleClass[iii] == i:
                        Number_Class += 1
                        if self.Align.RefList.loc[ii,iii] >0:
                            self.Align.RefList.at[ii,'Fill Group '+str(i)] += 1
                self.Align.RefList.loc[ii,'Fill Group '+str(i)] = self.Align.RefList.loc[ii,'Fill Group '+str(i)]/Number_Class
        self.Align.RefList['Fill QC'] = 0
        for ii in range(len(self.Align.RefList)):
            QC_array = self.Align.RefList.loc[ii,list(self.QCData.keys())]
            self.Align.RefList.at[ii,'Fill QC'] = len(QC_array[QC_array>=int(self.LineEdit_Int_min.text())])/len(list(self.QCData.keys()))
        #self.Align.RefList.to_excel(self.filepath_title+'Alignment-Origin-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        self.Align.BlankFilter()
        self.Align.RefList.insert(0,'ID',[0]*len(self.Align.RefList))
        self.Align.RefList.sort_values(by='m/z',ascending=True,inplace=True)
        self.Align.RefList.reset_index(drop=True,inplace=True)
        for ID_i in range(len(self.Align.RefList)):
            self.Align.RefList.at[ID_i,'ID'] = ID_i+1
            self.Align.RefList.at[ID_i,'m/z'] = np.around(self.Align.RefList['m/z'][ID_i],5)
            self.Align.RefList.at[ID_i,'RI'] = np.around(self.Align.RefList['RI'][ID_i],1)
            self.Align.RefList.at[ID_i,'RT'] = np.around(self.Align.RefList['RT'][ID_i],2)
        with open(self.filepath_title+'Align-RefList-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl','wb') as f:
            pickle.dump(self.Align.RefList,f)
        ''' Filter '''
        Filter_List = list(filter(lambda x:(self.Align.RefList.at[x,'Fill Group '+self.ClassList[0]]>=0.8 or self.Align.RefList.at[x,'Fill Group '+self.ClassList[1]]>=0.8) and (self.Align.RefList.at[x,'Fill QC']>=0.8),range(len(self.Align.RefList))))
        self.Align.RefList = self.Align.RefList.loc[Filter_List,:]
        self.Align.RefList.reset_index(drop=True,inplace=True)
        ''' MS2 assign'''
        for i in self.MS2_Data.keys():  
            V_List_Num = list(self.MS2_Data[i]['MS2_Data']['Pre_MZ'][int(len(self.MS2_Data[i]['MS2_Data']['Pre_MZ'])/2):int(len(self.MS2_Data[i]['MS2_Data']['Pre_MZ'])/2)+50])
            V_List_Count = []
            Count = 1
            for ii in range(1,len(V_List_Num)):
                if abs(V_List_Num[ii-1]-V_List_Num[ii])/V_List_Num[ii] <= float(self.LineEdit_MS_Tor.text())/1000000:
                    Count = Count + 1
                else:
                    V_List_Count.append(Count)
                    Count = 1
            V_List_Count = pd.Series(V_List_Count[2:len(V_List_Count)-2])
            V_List = []
            for ii in range(V_List_Count.mode()[0]):
                V_List.append('CE'+str(ii+1))
            break
        for i in self.MS2_Data.keys():
            self.MS2_Data[i]['MS2_Data']['CE'] = V_List[0]
            for ii in range(1,len(self.MS2_Data[i]['MS2_Data'])):
                if self.MS2_Data[i]['MS2_Data']['Pre_MZ'][ii] == self.MS2_Data[i]['MS2_Data']['Pre_MZ'][ii-1]:
                    V_Place = list(filter(lambda x:V_List[x]==self.MS2_Data[i]['MS2_Data']['CE'][ii-1],range(len(V_List))))[0]
                    if V_Place == len(V_List)-1:
                        self.MS2_Data[i]['MS2_Data'].loc[ii,'CE'] = V_List[0]
                    else:
                        self.MS2_Data[i]['MS2_Data'].loc[ii,'CE'] = V_List[V_Place+1]
                else:
                    self.MS2_Data[i]['MS2_Data'].loc[ii,'CE'] = V_List[0]
        self.MergeMS2 = []
        for i in self.MS2_Data.keys():
            if len(self.MergeMS2) == 0:
                self.MergeMS2 = self.MS2_Data[i]['MS2_Data'].copy()
            else:
                self.MergeMS2 = pd.concat([self.MergeMS2,self.MS2_Data[i]['MS2_Data']])
        if len(self.MergeMS2)>0:
            self.MergeMS2.reset_index(drop=True,inplace=True)
        self.Align.assign_MS2_KDtree(self.MergeMS2,float(self.LineEdit_RI_Tor.text())/100,float(self.LineEdit_MS_Tor.text())/1000000,V_List,RI_swich=True)
        ''' LOESS '''
        try:
            if len(list(self.QCData.keys()))>0:
                corrected_df = self.Align.RefList.copy()
                y_all = list(self.SampleData.keys())
                x_all = []
                for i in self.SampleData.keys():
                    x_all.append(float(self.SampleData[i]['Index']))
                x_all = np.array(x_all).reshape(-1, 1)
                qc_mask = list(self.QCData.keys())
                x_qc = []
                for i in self.QCData.keys():
                    x_qc.append(float(self.QCData[i]['Index']))
                x_qc = np.array(x_qc).reshape(-1, 1)
                for i in range(len(corrected_df)):
                    y_qc = list(corrected_df.loc[i,qc_mask])
                    y_qc = np.array(y_qc).reshape(-1, 1)
                    # 忽略常量/缺失
                    if np.std(y_qc) == 0 or len(y_qc[y_qc==0])/len(y_qc) > 0.2:
                        continue
                    else:
                        del_list = np.where(y_qc==0)
                        y_qc_fix = np.delete(y_qc,del_list[0])
                        y_qc_fix = np.array(y_qc_fix).reshape(-1, 1)
                        x_qc_fix = np.delete(x_qc,del_list[0])
                        x_qc_fix = np.array(x_qc_fix).reshape(-1, 1)
                    # 拟合 LASSO
                    loess_fit = lowess(endog=y_qc_fix.ravel(),
                       exog=x_qc_fix.ravel(),
                       frac=0.4, return_sorted=True)
                    x_fit, y_fit = loess_fit[:,0], loess_fit[:,1]
                    y_pred = np.interp(x_all.ravel(), x_fit, y_fit)
                    #model = ElasticNetCV(l1_ratio=0.5, cv=5).fit(x_qc_fix, y_qc_fix.ravel())
                    #y_pred = model.predict(x_all)
                    mean_qc = np.mean(y_qc)
                    corrected = list(corrected_df.loc[i,y_all] / y_pred * mean_qc)
                    corrected_df.loc[i,y_all] = corrected
                self.Origin_RefList = self.Align.RefList.copy()
                self.Align.RefList = corrected_df
                self.Align.RefList.loc[:,self.SampleData.keys()] = self.Align.RefList.loc[:,self.SampleData.keys()].astype(int)
        except Exception as e:
            print("Lasso回归异常",e)
        ''' --- Annotation --- '''
        try:
            def query_adducts_batch(db_path, table, mz_list, ri_list, mz_tol=0.000010, ri_tol=0.02):
                # table="adducts_positive" or "adducts_negative"
                conn = sqlite3.connect(db_path)

                conditions = []
                params = []  # 参数值
                ranges = []

                for mz, ri in zip(mz_list, ri_list):
                    mz_diff = mz * mz_tol
                    ri_diff = ri * ri_tol
                    mz_min, mz_max = mz - mz_diff, mz + mz_diff
                    ri_min, ri_max = ri - ri_diff, ri + ri_diff

                    conditions.append("(a.mz BETWEEN ? AND ? AND a.RI BETWEEN ? AND ?)")
                    params.extend([mz_min, mz_max, ri_min, ri_max])
                    ranges.append((mz_min, mz_max, ri_min, ri_max))

                query = f"""
                SELECT a.*, m."Compound Name", m."CAS Number", m."Compound ID", m."Formula", m."SMILES", m."InChIKey", m."KEGG"
                FROM {table} a
                JOIN metabolites m ON a.metabolite_id = m.id
                WHERE {" OR ".join(conditions)}
                """
                df_all = pd.read_sql_query(query, conn, params=params)
                conn.close()
                df_all = df_all.drop_duplicates(subset=["id"])
                df_all.reset_index(drop=True,inplace=True)
                df_all['MS2_MZ'] = df_all['MS2_MZ'].astype('object')
                df_all['MS2_Int'] = df_all['MS2_Int'].astype('object')

                for index in df_all.index:
                    if df_all.loc[index,'MS2_Int'] == 'NaN':
                        df_all.at[index,'MS2_MZ'] = []
                        df_all.at[index,'MS2_Int'] = []
                    else:
                        MS2_MZ = [float(x) for x in df_all.loc[index,'MS2_MZ'][2:-2].split(',')]
                        MS2_Int = [float(x) for x in df_all.loc[index,'MS2_Int'][2:-2].split(',')]
                        df_all.at[index,'MS2_MZ'] = MS2_MZ
                        df_all.at[index,'MS2_Int'] = MS2_Int
                # Python 拆分
                results = []
                for mz_min, mz_max, ri_min, ri_max in ranges:
                    mask = (
                        (df_all["mz"] >= mz_min) & (df_all["mz"] <= mz_max) &
                        (df_all["RI"] >= ri_min) & (df_all["RI"] <= ri_max)
                    )
                    results.append(df_all[mask].copy() if not df_all[mask].empty else pd.DataFrame())

                return results
            
            def query_adducts(db_path, table, mz_list, ri_list, Label_process_sub,process_bar, mz_tol=0.000010, ri_tol=0.02):
                # table="adducts_positive" or "adducts_negative"
                conn = sqlite3.connect(db_path)
                df_all = []
                query = f"""
                SELECT a.*, m."Compound Name", m."CAS Number", m."Compound ID", m."Formula", m."SMILES", m."InChIKey", m."KEGG"
                FROM {table} a
                JOIN metabolites m ON a.metabolite_id = m.id
                WHERE a.mz BETWEEN ? AND ? AND a.RI BETWEEN ? AND ?
                """
                Label_process_sub.setText('Annotation')
                process_bar.setMaximum(len(mz_list))
                process_bar.setValue(0)
                for mz, ri in zip(mz_list, ri_list):
                    process_bar.setValue(process_bar.value()+1)
                    QApplication.processEvents()
                    mz_diff = mz * mz_tol
                    ri_diff = ri * ri_tol
                    mz_min, mz_max = mz - mz_diff, mz + mz_diff
                    ri_min, ri_max = ri - ri_diff, ri + ri_diff
                    df = pd.read_sql_query(query, conn, params=(mz_min, mz_max,ri_min, ri_max))
                    df_all.append(df)
                conn.close()
                for i in range(len(df_all)):
                    if len(df_all[i]) > 0:
                        for ii in df_all[i].index:
                            df_all[i]['MS2_MZ'] = df_all[i]['MS2_MZ'].astype('object')
                            df_all[i]['MS2_Int'] = df_all[i]['MS2_Int'].astype('object')
                            if df_all[i].loc[ii,'MS2_Int'] == 'NaN':
                                df_all[i].at[ii,'MS2_MZ'] = []
                                df_all[i].at[ii,'MS2_Int'] = []
                            else:
                                MS2_MZ = [float(x) for x in df_all[i].loc[ii,'MS2_MZ'][2:-2].split(',')]
                                MS2_Int = [float(x) for x in df_all[i].loc[ii,'MS2_Int'][2:-2].split(',')]
                                df_all[i].at[ii,'MS2_MZ'] = MS2_MZ
                                df_all[i].at[ii,'MS2_Int'] = MS2_Int
                return df_all

            def MS2_Similarity(MS2_MZ_1, MS2_Int_1, MS2_MZ_2, MS2_Int_2, Number=5, MZ_Tor=0.000010):
                """
                优化版的MS2谱图相似性计算函数
                使用向量化操作和更高效的算法
                """
                # 转换为numpy数组以提高性能
                MS2_MZ_1 = np.asarray(MS2_MZ_1, dtype=np.float64)
                MS2_Int_1 = np.asarray(MS2_Int_1, dtype=np.float64)
                MS2_MZ_2 = np.asarray(MS2_MZ_2, dtype=np.float64)
                MS2_Int_2 = np.asarray(MS2_Int_2, dtype=np.float64)
                
                # 获取两个谱图中强度最高的Number个峰
                if len(MS2_Int_1) > Number:
                    top1_indices = np.argpartition(MS2_Int_1, -Number)[-Number:]
                    MZ_Standard_List = MS2_MZ_1[top1_indices]
                    Final_Int_1 = MS2_Int_1[top1_indices]
                else:
                    MZ_Standard_List = MS2_MZ_1.copy()
                    Final_Int_1 = MS2_Int_1.copy()
                
                if len(MS2_Int_2) > Number:
                    top2_indices = np.argpartition(MS2_Int_2, -Number)[-Number:]
                    MZ_append_List = MS2_MZ_2[top2_indices]
                    Int_append_List = MS2_Int_2[top2_indices]
                else:
                    MZ_append_List = MS2_MZ_2.copy()
                    Int_append_List = MS2_Int_2.copy()
                
                # 初始化第二个谱图的强度向量
                Final_Int_2 = np.zeros_like(Final_Int_1)
                
                # 第一步：为标准列表中的每个m/z在第二个谱图中寻找匹配
                for i, mz_std in enumerate(MZ_Standard_List):
                    # 计算相对质量偏差
                    rel_errors = np.abs(MS2_MZ_2 - mz_std) / mz_std
                    min_error_idx = np.argmin(rel_errors)
                    
                    if rel_errors[min_error_idx] < MZ_Tor:
                        Final_Int_2[i] = MS2_Int_2[min_error_idx]
                
                # 第二步：为第二个谱图中的高强峰在标准列表中寻找匹配或添加新峰
                new_mz_list = []
                new_int1_list = []
                new_int2_list = []
                
                for mz_app, int_app in zip(MZ_append_List, Int_append_List):
                    # 检查是否已经在标准列表中有匹配
                    rel_errors_std = np.abs(MZ_Standard_List - mz_app) / mz_app
                    min_std_error_idx = np.argmin(rel_errors_std)
                    
                    if rel_errors_std[min_std_error_idx] < MZ_Tor:
                        # 如果匹配到标准列表，更新强度（取较大值）
                        Final_Int_2[min_std_error_idx] = max(Final_Int_2[min_std_error_idx], int_app)
                    else:
                        # 检查是否在第一个谱图中有匹配
                        rel_errors_ms1 = np.abs(MS2_MZ_1 - mz_app) / mz_app
                        min_ms1_error_idx = np.argmin(rel_errors_ms1)
                        
                        if rel_errors_ms1[min_ms1_error_idx] < MZ_Tor:
                            # 添加到新峰列表
                            new_mz_list.append(mz_app)
                            new_int1_list.append(MS2_Int_1[min_ms1_error_idx])
                            new_int2_list.append(int_app)
                        else:
                            # 无匹配，第一个谱图强度为0
                            new_mz_list.append(mz_app)
                            new_int1_list.append(0.0)
                            new_int2_list.append(int_app)
                
                # 合并所有峰
                if new_mz_list:
                    Final_Int_1 = np.concatenate([Final_Int_1, new_int1_list])
                    Final_Int_2 = np.concatenate([Final_Int_2, new_int2_list])
                
                # 计算余弦相似度
                norm1 = np.linalg.norm(Final_Int_1)
                norm2 = np.linalg.norm(Final_Int_2)
                
                if norm1 == 0 or norm2 == 0:
                    return 0.0
                
                similarity = np.dot(Final_Int_1, Final_Int_2) / (norm1 * norm2)
                return similarity
            
            def export_network_to_graphml(RefList, Plt_DF, output_path='./Molecular_Network.graphml'):
                """
                将 MS/MS 相似性分子网络导出为 Cytoscape 可导入的 GraphML 文件
                ---------------------------------
                参数：
                    RefList : pd.DataFrame
                        包含 'm/z', 'RT', 'RI', 'Cor_Link', 'Cor_Value' 等列的 DataFrame
                    Plt_DF : pd.DataFrame
                        包含每个 networkx.Graph 的 DataFrame
                    output_path : str
                        输出文件路径（.graphml）
                """
                G_all = nx.Graph()
            
                for idx, row in Plt_DF.iterrows():
                    G = row['Graph']
                    if type(G) != list:
                        if len(G.nodes) == 0:
                            continue
                        network_id = f"Network_{idx}"
                
                        # 添加节点
                        for node in G.nodes():
                            mz = RefList.at[node, 'm/z']
                            rt = RefList.at[node, 'RT'] if 'RT' in RefList.columns else None
                            ri = RefList.at[node, 'RI'] if 'RI' in RefList.columns else None
                            Annotation = RefList.at[node, 'Annotation'] if 'Annotation' in RefList.columns else None
                            name = RefList.at[node, 'Name'] if 'Name' in RefList.columns else f"{node}"
                            if Annotation not in ["", "None", None]:
                                color = "#F29979"   # Annotated - Orange
                                category = "Annotated"
                            else:
                                color = "#C0DEED"   # Unknow - Blue B0BEC5
                                category = "Unannotated"
                            G_all.add_node(node,
                                Name=name,
                                mz=float(mz),
                                RT=float(rt) if rt is not None else None,
                                RI=float(ri) if ri is not None else None,
                                Annotation = Annotation,
                                Network=network_id,
                                Color=color,           # Cytoscape 可识别的节点颜色
                                Category=category      # 方便分组筛选
                            )
                
                        # 添加边
                        for u, v, d in G.edges(data=True):
                            weight = d.get('weight', 0)
                            G_all.add_edge(u, v, weight=float(weight), Network=network_id)
            
                # 输出 GraphML
                nx.write_graphml(G_all, output_path)
                return output_path

            if self.Polarity == 'Positive':
                table = 'adducts_positive'
            else:
                table = 'adducts_negative'
            mz_list = list(self.Align.RefList['m/z'])
            ri_list = list(self.Align.RefList['RI'])
            Annotation = query_adducts(self.DB_Path, table, mz_list, ri_list, self.Label_process_sub, self.process_bar, mz_tol=float(self.LineEdit_MS_Tor.text())/1000000, ri_tol=float(self.LineEdit_RI_Tor.text())/100)
            self.Align.RefList['Annotation'] = ''
            for i in range(len(Annotation)):
                if len(Annotation[i]) == 0:
                    continue
                else:
                    Annotate_results = ''
                    if len(self.Align.RefList.loc[i,'MS2_MZ'])>0:
                        for ii in Annotation[i].index:
                            if len(Annotation[i].loc[ii,'MS2_MZ'])>0:
                                match_score = MS2_Similarity(self.Align.RefList.loc[i,'MS2_MZ'],self.Align.RefList.loc[i,'MS2_Int'],Annotation[i].loc[ii,'MS2_MZ'],Annotation[i].loc[ii,'MS2_Int'],Number=int(self.LineEdit_Annotation_MS2_Number.text()),MZ_Tor=float(self.LineEdit_MS_Tor.text())/1000000)
                                if match_score > float(self.LineEdit_min_cos.text()):
                                    if len(Annotate_results) == 0:
                                        Annotate_results += '[m/z+RI+Tandem MS:'+'Similarity Score:'+str(np.around(match_score,3))+';Name:'+Annotation[i].loc[ii,'Compound Name']+';CAS:'+str(Annotation[i].loc[ii,'CAS Number'])+';Compound ID:'+str(Annotation[i].loc[ii,'Compound ID'])+';Formula:'+str(Annotation[i].loc[ii,'Formula'])+';SMILES:'+str(Annotation[i].loc[ii,'SMILES'])+';InChIKey:'+str(Annotation[i].loc[ii,'InChIKey'])+';KEGG:'+str(Annotation[i].loc[ii,'KEGG'])+']'
                                    else:
                                        Annotate_results += ',[m/z+RI+Tandem MS:'+'Similarity Score:'+str(np.around(match_score,3))+';Name:'+Annotation[i].loc[ii,'Compound Name']+';CAS:'+str(Annotation[i].loc[ii,'CAS Number'])+';Compound ID:'+str(Annotation[i].loc[ii,'Compound ID'])+';Formula:'+str(Annotation[i].loc[ii,'Formula'])+';SMILES:'+str(Annotation[i].loc[ii,'SMILES'])+';InChIKey:'+str(Annotation[i].loc[ii,'InChIKey'])+';KEGG:'+str(Annotation[i].loc[ii,'KEGG'])+']'
                            else:
                                if len(Annotate_results) == 0:
                                    Annotate_results += '[m/z+RI: Name:'+Annotation[i].loc[ii,'Compound Name']+';CAS:'+str(Annotation[i].loc[ii,'CAS Number'])+';Compound ID:'+str(Annotation[i].loc[ii,'Compound ID'])+';Formula:'+str(Annotation[i].loc[ii,'Formula'])+';SMILES:'+str(Annotation[i].loc[ii,'SMILES'])+';InChIKey:'+str(Annotation[i].loc[ii,'InChIKey'])+';KEGG:'+str(Annotation[i].loc[ii,'KEGG'])+']'
                                else:
                                    Annotate_results += ',[m/z+RI: Name:'+Annotation[i].loc[ii,'Compound Name']+';CAS:'+str(Annotation[i].loc[ii,'CAS Number'])+';Compound ID:'+str(Annotation[i].loc[ii,'Compound ID'])+';Formula:'+str(Annotation[i].loc[ii,'Formula'])+';SMILES:'+str(Annotation[i].loc[ii,'SMILES'])+';InChIKey:'+str(Annotation[i].loc[ii,'InChIKey'])+';KEGG:'+str(Annotation[i].loc[ii,'KEGG'])+']'
                    else:
                        for ii in Annotation[i].index:
                            if len(Annotate_results) == 0:
                                Annotate_results += '[m/z+RI: Name:'+Annotation[i].loc[ii,'Compound Name']+';CAS:'+str(Annotation[i].loc[ii,'CAS Number'])+';Compound ID:'+str(Annotation[i].loc[ii,'Compound ID'])+';Formula:'+str(Annotation[i].loc[ii,'Formula'])+';SMILES:'+str(Annotation[i].loc[ii,'SMILES'])+';InChIKey:'+str(Annotation[i].loc[ii,'InChIKey'])+';KEGG:'+str(Annotation[i].loc[ii,'KEGG'])+']'
                            else:
                                Annotate_results += ',[m/z+RI: Name:'+Annotation[i].loc[ii,'Compound Name']+';CAS:'+str(Annotation[i].loc[ii,'CAS Number'])+';Compound ID:'+str(Annotation[i].loc[ii,'Compound ID'])+';Formula:'+str(Annotation[i].loc[ii,'Formula'])+';SMILES:'+str(Annotation[i].loc[ii,'SMILES'])+';InChIKey:'+str(Annotation[i].loc[ii,'InChIKey'])+';KEGG:'+str(Annotation[i].loc[ii,'KEGG'])+']'
                self.Align.RefList.at[i,'Annotation'] = Annotate_results
            with open(self.filepath_title+'Align-RefList-Annotation-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl','wb') as f:
                pickle.dump(self.Align.RefList,f)
            Annotate_results = []
            for i in range(len(Annotation)):
                if len(Annotation[i]) == 0:
                    continue
                else:
                    for ii in Annotation[i].index:
                        if len(Annotation[i].loc[ii,'MS2_MZ'])>0 and len(self.Align.RefList.loc[i,'MS2_MZ'])>0:
                            match_score = MS2_Similarity(self.Align.RefList.loc[i,'MS2_MZ'],self.Align.RefList.loc[i,'MS2_Int'],Annotation[i].loc[ii,'MS2_MZ'],Annotation[i].loc[ii,'MS2_Int'],Number=5,MZ_Tor=float(self.LineEdit_MS_Tor.text())/1000000)
                            if match_score > 0.7:
                                temp_df = self.Align.RefList.loc[[i],:]
                                temp_df['Match Type'] = 'm/z+RI+Tandem MS'
                                temp_df['Similarity Score'] = np.around(match_score,3)
                                temp_df['Name'] = Annotation[i].loc[ii,'Compound Name']
                                temp_df['CAS Number'] = str(Annotation[i].loc[ii,'CAS Number'])
                                temp_df['Compound ID'] = str(Annotation[i].loc[ii,'Compound ID'])
                                temp_df['Formula'] = str(Annotation[i].loc[ii,'Formula'])
                                temp_df['SMILES'] = str(Annotation[i].loc[ii,'SMILES'])
                                temp_df['InChIKey'] = str(Annotation[i].loc[ii,'InChIKey'])
                                temp_df['KEGG'] = str(Annotation[i].loc[ii,'KEGG'])
                                Annotate_results.append(temp_df)
                        else:
                            temp_df = self.Align.RefList.loc[[i],:]
                            temp_df['Match Type'] = 'm/z+RI'
                            temp_df['Similarity Score'] = ''
                            temp_df['Name'] = Annotation[i].loc[ii,'Compound Name']
                            temp_df['CAS Number'] = str(Annotation[i].loc[ii,'CAS Number'])
                            temp_df['Compound ID'] = str(Annotation[i].loc[ii,'Compound ID'])
                            temp_df['Formula'] = str(Annotation[i].loc[ii,'Formula'])
                            temp_df['SMILES'] = str(Annotation[i].loc[ii,'SMILES'])
                            temp_df['InChIKey'] = str(Annotation[i].loc[ii,'InChIKey'])
                            temp_df['KEGG'] = str(Annotation[i].loc[ii,'KEGG'])
                            Annotate_results.append(temp_df)
            self.Annotate_results = pd.concat(Annotate_results,ignore_index=True)
            self.Annotate_results.to_excel(self.filepath_title+'AnnotationResults-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
            '''--------- Network -------------'''
            self.Align.RefList['Cor_Link']=self.Align.RefList['m/z'].apply(lambda x:[])
            self.Align.RefList['Cor_Value']=self.Align.RefList['m/z'].apply(lambda x:[])
            self.Label_process_sub.setText('Network')
            self.process_bar.setMaximum(len(self.Align.RefList))
            self.process_bar.setValue(0)
            QApplication.processEvents()
            for i in range(len(self.Align.RefList)):
                self.process_bar.setValue(i)
                QApplication.processEvents()
                if len(self.Align.RefList.loc[i,'MS2_Int'])>0:
                    for ii in range(i+1,len(self.Align.RefList)):
                        if len(self.Align.RefList.loc[ii,'MS2_Int'])>0:
                            MSList1 = self.Align.RefList.at[i,'MS2_MZ']
                            IntList1 = self.Align.RefList.at[i,'MS2_Int']
                            MSList2 = self.Align.RefList.at[ii,'MS2_MZ']
                            IntList2 = self.Align.RefList.at[ii,'MS2_Int']
                            MSMS_Similarity = MS2_Similarity(MSList1,IntList1,MSList2,IntList2,Number=5,MZ_Tor=0.000010)
                            if MSMS_Similarity > 0.7:
                                self.Align.RefList.at[i,'Cor_Link'].append(ii)
                                self.Align.RefList.at[ii,'Cor_Link'].append(i)
                                self.Align.RefList.at[i,'Cor_Value'].append(MSMS_Similarity)
                                self.Align.RefList.at[ii,'Cor_Value'].append(MSMS_Similarity)
            Plt_index = list(filter(lambda x:len(self.Align.RefList.loc[x,'Cor_Link'])>0,range(len(self.Align.RefList))))
            Plt_DF = pd.DataFrame({'Rank':range(len(Plt_index))})
            Plt_DF['Graph']=Plt_DF['Rank'].apply(lambda x:[])
            bar = Bar('Draw Molecule Network', max=len(self.Align.RefList))
            for i in Plt_index:
                bar.next() 
                Exest = list(filter(lambda x:i in Plt_DF.loc[x,'Graph'],range(len(Plt_DF))))
                if len(Exest)==0:
                    G = nx.Graph()
                    for ii in range(len(self.Align.RefList.loc[i,'Cor_Link'])):
                        G.add_edge(i,self.Align.RefList.loc[i,'Cor_Link'][ii],weight=self.Align.RefList.loc[i,'Cor_Value'][ii])
                    To_add_G = list(filter(lambda x:len(Plt_DF['Graph'][x])==0,range(len(Plt_DF))))
                    Plt_DF.at[To_add_G[0],'Graph']=G
                elif len(Exest)==1:
                    G = Plt_DF.loc[Exest[0],'Graph']
                    for ii in range(len(self.Align.RefList.loc[i,'Cor_Link'])):
                        G.add_edge(i,self.Align.RefList.loc[i,'Cor_Link'][ii],weight=self.Align.RefList.loc[i,'Cor_Value'][ii])
                    Plt_DF.at[Exest[0],'Graph']=G
            export_network_to_graphml(self.Align.RefList,Plt_DF, output_path=self.filepath_title+'Molecular_Network.graphml')
        except Exception as e:
            print("Annotation error",e)
        self.RIIS = pd.read_excel(main_HP.RIISPath)
        self.RIIS['RTList'] = self.RIIS['C'].apply(lambda x:[])
        Index = []
        Name_RIIS = []
        for i in self.SampleData.keys():
            Index.append(float(self.SampleData[i]['Index']))
            Name_RIIS.append(self.SampleData[i]['Name'])
            for ii in range(len(self.SampleData[i]['RIIS'])):
                self.RIIS.at[ii,'RTList'].append(self.SampleData[i]['RIIS'].at[ii,'RT'])
                
        for i in self.QCData.keys():
            Index.append(float(self.QCData[i]['Index']))
            Name_RIIS.append(self.QCData[i]['Name'])
            for ii in range(len(self.QCData[i]['RIIS'])):
                self.RIIS.at[ii,'RTList'].append(self.QCData[i]['RIIS'].at[ii,'RT'])
                
        fig = go.Figure()
        for ii in range(len(self.RIIS)):
            fig.add_trace(
                go.Scatter(
                    x=Index,
                    y=np.array(self.RIIS.at[ii,'RTList'])/60,
                    customdata=Name_RIIS,  # 绑定额外数据
                    hovertemplate=(
                    "Index: %{x}<br>"
                    "RT: %{y}<br>"
                    "Sample: %{customdata}<br>"  # 引用 customdata
                    "<extra></extra>"),
                    mode='markers',
                    name='C-'+str(self.RIIS.at[ii,'C']),
                    line={'width':2},
                    ))
        fig.update_layout(
            legend={
                'xanchor': 'right',
                'x': 1.1,  # 改为 0.95，让图例在画布内部（95% 宽度处）
                'y': 0.5,   # 垂直居中
                'bgcolor': 'rgba(255,255,255,0.7)',  # 可选：添加半透明背景
            },
            margin=dict(l=50, r=150, b=50, t=50, pad=4),
            autosize=False,
            width=1000,
            height=700,
            title='Index - RT of calibrants',
            titlefont={'size':20},
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis={'title':{'text':'Index','font':{'size':15},'standoff':0},
                   'linecolor':'black',
                   'tickfont':{'size':11},
                   'ticks':'outside',
                   'ticklen':2,
                   # 添加range设置，在两端留出空白
                   #'range':[-0.5, len(x_List)-0.5]  # 这里-0.5和+0.5表示两端各留出相当于半个类别的空白
                   },
            yaxis={'title':{'text':'RT(min)','font':{'size':15},'standoff':0},
                   'linecolor':'black',
                   'tickfont':{'size':11},
                   'ticks':'outside',
                   'ticklen':2,
                   'side':'left',
                   #'range':(0,1500),
                   }
            )
        fig.layout.font.family = 'Helvetica'
        fig.write_html(self.filepath_title+'Calibrants-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.html',config={'responsive': False})
        ''' --- statistics ---'''    
        def calculate_vip(pls):
            t = pls.x_scores_
            w = pls.x_weights_
            q = pls.y_loadings_
            p,h = w.shape
            s = np.sum(t ** 2,axis=0) * (q[0] ** 2)
            total_s = np.sum(s)
            vip = np.zeros(p)
            for j in range(p):
                weight_sum = 0
                for i in range(h):
                    weight_sum += s[i] * (w[j,i]**2)
                vip[j] = np.sqrt(p * weight_sum / total_s)
            return vip
        try:
            Filter_List = list(filter(lambda x:(self.Align.RefList.at[x,'Fill Group '+self.ClassList[0]]>=0.8 or self.Align.RefList.at[x,'Fill Group '+self.ClassList[1]]>=0.8) and (self.Align.RefList.at[x,'Fill QC']>=0.8),range(len(self.Align.RefList))))
            self.Filter_RefList = self.Align.RefList.loc[Filter_List,:]
            self.Filter_RefList.reset_index(drop=True,inplace=True)
            Filter_RefList_Fix = self.Filter_RefList
            Group_1_samples = [s for s in list(self.SampleData.keys()) if self.SampleData[s]['Group'] == self.ClassList[0]]
            Group_1_samples.sort()
            Group_2_samples = [s for s in list(self.SampleData.keys()) if self.SampleData[s]['Group'] == self.ClassList[1]]
            Group_2_samples.sort()
            self.Class_dict = {}
            Class_All = []
            for i in self.ClassList:
                self.Class_dict[i] = [s for s in list(self.SampleData.keys()) if self.SampleData[s]['Group'] == i]
                self.Class_dict[i].sort()
                Class_All += self.Class_dict[i]
            Order = [x for x in range(len(list(self.Filter_RefList.keys()))) if list(self.Filter_RefList.keys())[x] in Class_All]
            Key_List = np.array(self.Filter_RefList.keys())
            Key_List[Order] = Class_All
            self.Filter_RefList = self.Filter_RefList.loc[:,Key_List]
            QC_samples = list(self.QCData.keys())
            # Lod and max values
            for i in range(len(Filter_RefList_Fix)):
                for ii in list(self.Class_dict.keys()):
                    p95 = Filter_RefList_Fix.loc[i,self.Class_dict[ii]].quantile(0.95)
                    Filter_RefList_Fix.loc[i,self.Class_dict[ii]] = Filter_RefList_Fix.loc[i,self.Class_dict[ii]].clip(lower=int(self.LineEdit_Int_min.text())/5,upper=p95)
            p95 = Filter_RefList_Fix.loc[i,QC_samples].quantile(0.95)
            Filter_RefList_Fix.loc[i,QC_samples] = Filter_RefList_Fix.loc[i,QC_samples].clip(lower=int(self.LineEdit_Int_min.text())/5,upper=p95)
            # PCA
            if len(self.ClassList) == 2: 
                X = Filter_RefList_Fix.loc[:,Group_1_samples+Group_2_samples+QC_samples].T
                X_scaled = StandardScaler().fit_transform(X)
                pca = PCA(n_components=2)
                scores = pca.fit_transform(X_scaled)
                explained_var = pca.explained_variance_ratio_ * 100  # 转换为百分比
                fig = go.Figure()
                # Group 1
                data = scores[range(len(Group_1_samples)), :2]
                fig.add_trace(
                    go.Scatter(
                        x=data[:, 0],
                        y=data[:, 1],
                        mode='markers',
                        name=self.ClassList[0],
                        line={'width':2},
                        marker=dict(
                            size=10,
                            color='rgba(247,174,177,1)',
                            line=dict(width=1, color="black")
                        )))
                mean = np.mean(data, axis=0)
                cov = np.cov(data, rowvar=False)
                chi2_val = chi2.ppf(0.95, df=2)
            
                eigvals, eigvecs = np.linalg.eigh(cov)
                order = eigvals.argsort()[::-1]
                eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
                theta = np.linspace(0, 2*np.pi, 200)
                ellipse = np.array([np.cos(theta), np.sin(theta)])
                ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                ellipse = ellipse.T + mean
            
                fig.add_trace(go.Scatter(
                    x=ellipse[:, 0],
                    y=ellipse[:, 1],
                    mode="lines",
                    line=dict(color='rgba(247,174,177,1)', width=2),
                    fill="toself",
                    fillcolor='rgba(247,174,177,0.3)',
                    name=f"{self.ClassList[0]} 95% CI"
                ))
                # Group 2
                data = scores[range(len(Group_1_samples),len(Group_1_samples)+len(Group_2_samples)), :2]
                fig.add_trace(
                    go.Scatter(
                        x=data[:, 0],
                        y=data[:, 1],
                        mode='markers',
                        name=self.ClassList[1],
                        line={'width':2},
                        marker=dict(
                            size=10,
                            color='rgba(174,204,231,1)',
                            line=dict(width=1, color="black")
                        )))
                mean = np.mean(data, axis=0)
                cov = np.cov(data, rowvar=False)
                chi2_val = chi2.ppf(0.95, df=2)
            
                eigvals, eigvecs = np.linalg.eigh(cov)
                order = eigvals.argsort()[::-1]
                eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
                theta = np.linspace(0, 2*np.pi, 200)
                ellipse = np.array([np.cos(theta), np.sin(theta)])
                ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                ellipse = ellipse.T + mean
            
                fig.add_trace(go.Scatter(
                    x=ellipse[:, 0],
                    y=ellipse[:, 1],
                    mode="lines",
                    line=dict(color='rgba(174,204,231,1)', width=2),
                    fill="toself",
                    fillcolor='rgba(174,204,231,0.3)',
                    name=f"{self.ClassList[1]} 95% CI"
                ))
                # Group QC
                data = scores[range(len(Group_1_samples)+len(Group_2_samples),len(Group_1_samples)+len(Group_2_samples)+len(QC_samples)), :2]
                fig.add_trace(
                    go.Scatter(
                        x=data[:, 0],
                        y=data[:, 1],
                        mode='markers',
                        name='QC',
                        line={'width':2},
                        marker=dict(
                            size=10,
                            color='rgba(143,180,137,1)',
                            line=dict(width=1, color="black")
                        )))
                mean = np.mean(data, axis=0)
                cov = np.cov(data, rowvar=False)
                chi2_val = chi2.ppf(0.95, df=2)
            
                eigvals, eigvecs = np.linalg.eigh(cov)
                order = eigvals.argsort()[::-1]
                eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
                theta = np.linspace(0, 2*np.pi, 200)
                ellipse = np.array([np.cos(theta), np.sin(theta)])
                ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                ellipse = ellipse.T + mean
            
                fig.add_trace(go.Scatter(
                    x=ellipse[:, 0],
                    y=ellipse[:, 1],
                    mode="lines",
                    line=dict(color='rgba(143,180,137,1)', width=2),
                    fill="toself",
                    fillcolor='rgba(143,180,137,0.3)',
                    name=f"{'QC'} 95% CI"
                ))
                fig.update_layout(
                    width=700,
                    height=600,
                    title="PCA 2D Scores Plot",
                    xaxis_title=f"PC1 ({explained_var[0]:.2f}%)",
                    yaxis_title=f"PC2 ({explained_var[1]:.2f}%)",
                    template="simple_white",
                    legend=dict(itemsizing="constant")
                    )
                fig.write_html(self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.html',config={'responsive': False})
                #pio.write_image(fig,self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.svg')
                data = pd.DataFrame({'Sample':Group_1_samples+Group_2_samples+QC_samples,'PC1':scores[:,0],'PC2':scores[:,1]})
                data.to_excel(self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
                with open(self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.txt','w') as f:
                    f.write('PC1 explained_variance_ratio:'+str(explained_var[0])+'\nPC2 explained_variance_ratio:'+str(explained_var[1]))
            else:
                X = Filter_RefList_Fix.loc[:,Class_All+QC_samples].T
                X_scaled = StandardScaler().fit_transform(X)
                pca = PCA(n_components=2)
                scores = pca.fit_transform(X_scaled)
                explained_var = pca.explained_variance_ratio_ * 100  # 转换为百分比
                range_start = 0
                colors = px.colors.qualitative.Plotly
                fig = go.Figure()
                for i in range(len(self.ClassList)):
                    base_color = colors[i % len(colors)]
                    hex_color = base_color.lstrip('#')
                    rgb = tuple(int(hex_color[i:i+2],16) for i in (0,2,4))
                    line_color = f'rgb{rgb}'
                    fill_color = f'rgba({rgb[0]},{rgb[1]},{rgb[2]},0.2)'
                    data = scores[range(range_start,range_start+len(self.Class_dict[self.ClassList[i]])), :2]
                    fig.add_trace(
                        go.Scatter(
                            x=data[:, 0],
                            y=data[:, 1],
                            mode='markers',
                            name=self.ClassList[i],
                            line={'width':2},
                            marker=dict(
                                size=10,
                                color=line_color,
                                line=dict(width=1, color="black")
                            )))
                    mean = np.mean(data, axis=0)
                    cov = np.cov(data, rowvar=False)
                    chi2_val = chi2.ppf(0.95, df=2)
                
                    eigvals, eigvecs = np.linalg.eigh(cov)
                    order = eigvals.argsort()[::-1]
                    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
                
                    theta = np.linspace(0, 2*np.pi, 200)
                    ellipse = np.array([np.cos(theta), np.sin(theta)])
                    ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                    ellipse = ellipse.T + mean
                    
                    fig.add_trace(go.Scatter(
                        x=ellipse[:, 0],
                        y=ellipse[:, 1],
                        mode="lines",
                        line=dict(color=line_color, width=2),
                        fill="toself",
                        fillcolor=fill_color,
                        name=f"{self.ClassList[i]} 95% CI"
                    ))
                    range_start += len(self.Class_dict[self.ClassList[i]])
                # Group QC
                base_color = colors[i+1 % len(colors)]
                hex_color = base_color.lstrip('#')
                rgb = tuple(int(hex_color[i:i+2],16) for i in (0,2,4))
                line_color = f'rgb{rgb}'
                fill_color = f'rgba({rgb[0]},{rgb[1]},{rgb[2]},0.2)'
                data = scores[range(range_start,range_start+len(QC_samples)), :2]
                fig.add_trace(
                    go.Scatter(
                        x=data[:, 0],
                        y=data[:, 1],
                        mode='markers',
                        name='QC',
                        line={'width':2},
                        marker=dict(
                            size=10,
                            color=line_color,
                            line=dict(width=1, color="black")
                        )))
                mean = np.mean(data, axis=0)
                cov = np.cov(data, rowvar=False)
                chi2_val = chi2.ppf(0.95, df=2)
            
                eigvals, eigvecs = np.linalg.eigh(cov)
                order = eigvals.argsort()[::-1]
                eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
                theta = np.linspace(0, 2*np.pi, 200)
                ellipse = np.array([np.cos(theta), np.sin(theta)])
                ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                ellipse = ellipse.T + mean
            
                fig.add_trace(go.Scatter(
                    x=ellipse[:, 0],
                    y=ellipse[:, 1],
                    mode="lines",
                    line=dict(color=line_color, width=2),
                    fill="toself",
                    fillcolor=fill_color,
                    name=f"{'QC'} 95% CI"
                ))
                fig.update_layout(
                    width=700,
                    height=600,
                    title="PCA 2D Scores Plot",
                    xaxis_title=f"PC1 ({explained_var[0]:.2f}%)",
                    yaxis_title=f"PC2 ({explained_var[1]:.2f}%)",
                    template="simple_white",
                    legend=dict(itemsizing="constant")
                    )
                fig.write_html(self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.html',config={'responsive': False})
                #pio.write_image(fig,self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.svg')
                data = pd.DataFrame({'Sample':Class_All,'PC1':scores[:,0],'PC2':scores[:,1]})
                data.to_excel(self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
                with open(self.filepath_title+'PCA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.txt','w') as f:
                    f.write('PC1 explained_variance_ratio:'+str(explained_var[0])+'\nPC2 explained_variance_ratio:'+str(explained_var[1]))
            # PLS-DA
            if len(self.ClassList) == 2:
                y_origin = [1]*len(Group_1_samples)+[0]*len(Group_2_samples)
                y_enc = LabelEncoder().fit_transform(y_origin)
                scaler = StandardScaler()
                x_origin = Filter_RefList_Fix.loc[:,Group_1_samples+Group_2_samples]
                x_origin = x_origin.T.values
                x_scaled = scaler.fit_transform(x_origin)
                pls = PLSRegression(n_components=2)
                pls.fit(x_scaled, y_enc)
                vip_scores = calculate_vip(pls)
                scores = pls.x_scores_
                self.Filter_RefList['PLS-DA-VIP'] = vip_scores
                Filter_RefList_Fix['PLS-DA-VIP'] = vip_scores
                fig = go.Figure()
                # Group 1
                data = scores[range(len(Group_1_samples)), :2]
                fig.add_trace(
                    go.Scatter(
                        x=data[:, 0],
                        y=data[:, 1],
                        mode='markers',
                        name=self.ClassList[0],
                        line={'width':2},
                        marker=dict(
                            size=10,
                            color='rgba(247,174,177,1)',
                            line=dict(width=1, color="black")
                        )))
                mean = np.mean(data, axis=0)
                cov = np.cov(data, rowvar=False)
                chi2_val = chi2.ppf(0.95, df=2)
            
                eigvals, eigvecs = np.linalg.eigh(cov)
                order = eigvals.argsort()[::-1]
                eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
                theta = np.linspace(0, 2*np.pi, 200)
                ellipse = np.array([np.cos(theta), np.sin(theta)])
                ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                ellipse = ellipse.T + mean
                
                fig.add_trace(go.Scatter(
                    x=ellipse[:, 0],
                    y=ellipse[:, 1],
                    mode="lines",
                    line=dict(color='rgba(247,174,177,1)', width=2),
                    fill="toself",
                    fillcolor='rgba(247,174,177,0.3)',
                    name=f"{self.ClassList[0]} 95% CI"
                ))
                # Group 2
                data = scores[range(len(Group_1_samples),len(Group_1_samples)+len(Group_2_samples)), :2]
                fig.add_trace(
                    go.Scatter(
                        x=data[:, 0],
                        y=data[:, 1],
                        mode='markers',
                        name=self.ClassList[1],
                        line={'width':2},
                        marker=dict(
                            size=10,
                            color='rgba(174,204,231,1)',
                            line=dict(width=1, color="black")
                        )))
                mean = np.mean(data, axis=0)
                cov = np.cov(data, rowvar=False)
                chi2_val = chi2.ppf(0.95, df=2)
            
                eigvals, eigvecs = np.linalg.eigh(cov)
                order = eigvals.argsort()[::-1]
                eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
                theta = np.linspace(0, 2*np.pi, 200)
                ellipse = np.array([np.cos(theta), np.sin(theta)])
                ellipse = np.sqrt(chi2_val) * eigvecs @ np.diag(np.sqrt(eigvals)) @ ellipse
                ellipse = ellipse.T + mean
            
                fig.add_trace(go.Scatter(
                    x=ellipse[:, 0],
                    y=ellipse[:, 1],
                    mode="lines",
                    line=dict(color='rgba(174,204,231,1)', width=2),
                    fill="toself",
                    fillcolor='rgba(174,204,231,0.3)',
                    name=f"{self.ClassList[1]} 95% CI"
                ))
                fig.update_layout(
                    width=700,
                    height=600,
                    title="PLS-DA 2D Scores Plot",
                    xaxis_title='LV1',
                    yaxis_title='LV2',
                    template="simple_white",
                    legend=dict(itemsizing="constant")
                    )
                fig.write_html(self.filepath_title+'PLS-DA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.html',config={'responsive': False})
                #pio.write_image(fig,self.filepath_title+'PLS-DA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.svg')
                data = pd.DataFrame({'Sample':Group_1_samples+Group_2_samples,'LV1':scores[:,0],'LV2':scores[:,1]})
                data.to_excel(self.filepath_title+'PLS-DA-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
            # Fold Change
            FC_list = []
            pvalue_list = []
            for peak in Filter_RefList_Fix.index:
                Values_1 = Filter_RefList_Fix.loc[peak,Group_1_samples].astype(float)
                Values_2 = Filter_RefList_Fix.loc[peak,Group_2_samples].astype(float)
                mean_1 = np.mean(Values_1)
                mean_2 = np.mean(Values_2)
                if mean_1 == 0 :
                    FC = np.nan
                else:
                    FC = mean_2/mean_1
                FC_list.append(FC)
                stat,pvalue = ttest_ind(Values_2,Values_1,equal_var=False)
                pvalue_list.append(pvalue)
            self.Filter_RefList['FC'] = FC_list
            Filter_RefList_Fix['FC'] = FC_list
            reject, p_adj, _, _ = multipletests(pvalue_list, method="fdr_bh")
            self.Filter_RefList['FC-p values'] = p_adj
            Filter_RefList_Fix['FC-p values'] = p_adj
            self.Filter_RefList.reset_index(drop=True,inplace=True)
            self.Align.RefList = self.Filter_RefList.copy()
            #self.Filter_RefList = self.Filter_RefList[((self.Filter_RefList['FC']>=2) | (self.Filter_RefList['FC']<=0.5)) & (self.Filter_RefList['FC-p values']<=0.05) & (self.Filter_RefList['PLS-DA-VIP']>=1)]
            self.Filter_RefList.reset_index(drop=True,inplace=True)
            
        except Exception as e:
            print("Statistics error",e)
        
        ''' --- Output Filtered Alignment List ---'''
        try:
            mgf_output = ['0']*1000000
            mgf_place = 0
            csv_mgf = pd.DataFrame(columns=['row ID','row m/z','row retention time','Peak height'])
            for i in  range(len(self.Filter_RefList)):
                if type(self.Filter_RefList['MS2_MZ'][i]) == list:
                    MZ = self.Filter_RefList.at[i,'m/z']
                    RT = np.around(self.Filter_RefList.at[i,'RT']/60,2)
                    Int = self.Filter_RefList['Int'][i]     
                    temp_csv_mgf = pd.DataFrame([(i+1,MZ,RT,Int)],columns=['row ID','row m/z','row retention time','Peak height'])
                    csv_mgf = pd.concat([csv_mgf,temp_csv_mgf])
                    if len(self.Filter_RefList.at[i,'MS2_Int']) > 0:
                        if self.Polarity == 'Positive':
                            if mgf_place >= len(mgf_output):
                                mgf_output = mgf_output + ['0']*1000000
                            mgf_output[mgf_place] = 'BEGIN IONS\nFEATURE_ID='+str(i+1)+'\nPEPMASS='+str(MZ)+'\nSCANS='+str(i+1)+'\nRTINSECONDS='+str(self.Filter_RefList.at[i,'RT'])+'\nCHARGE=+1\nMSLEVEL=2\n'
                            mgf_place += 1
                        elif self.Polarity == 'Negative':
                            if mgf_place >= len(mgf_output):
                                mgf_output = mgf_output + ['0']*1000000
                            mgf_output[mgf_place] = 'BEGIN IONS\nFEATURE_ID='+str(i+1)+'\nPEPMASS='+str(MZ)+'\nSCANS='+str(i+1)+'\nRTINSECONDS='+str(self.Filter_RefList.at[i,'RT'])+'\nCHARGE=-1\nMSLEVEL=2\n'
                            mgf_place += 1
                        for ii in range(len(self.Filter_RefList.at[i,'MS2_MZ'])):
                            if mgf_place >= len(mgf_output):
                                mgf_output = mgf_output + ['0']*1000000
                            mgf_output[mgf_place] = str(np.around(self.Filter_RefList.at[i,'MS2_MZ'][ii],4))+'\t'+str(int(self.Filter_RefList.at[i,'MS2_Int'][ii]))+'\n'
                            mgf_place += 1
                        if mgf_place >= len(mgf_output):
                            mgf_output = mgf_output + ['0']*1000000
                        mgf_output[mgf_place] = 'END IONS\n\n'
                        mgf_place += 1
            self.mgf_output = np.array(mgf_output)
            for i_mgf in range(len(mgf_output)-1,-1,-1):
                if mgf_output[i_mgf] != '0':
                    break
            mgf_output = mgf_output[0:i_mgf+1]
            mgf_str = ''.join(mgf_output)      
            if len(mgf_str)>0:
                with open(self.filepath_title+'tandemMS-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.mgf','w')as mgfFile:
                    mgfFile.write(mgf_str)
                csv_mgf.to_csv(self.filepath_title+'tandemMS-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.csv',index=False)
        except Exception as e:
            print("Output异常",e)
        Filter_RefList_output = self.Filter_RefList.copy()
        Filter_RefList_output.rename(columns={'Int':'Average Intensity'},inplace=True)
        Filter_RefList_output.rename(columns={'RT':'Average RT(s)'},inplace=True)
        Filter_RefList_output.rename(columns={'RI':'Average RI'},inplace=True)
        Filter_RefList_output.rename(columns={'MS_List':'m/z value in each sample'},inplace=True)
        Filter_RefList_output.rename(columns={'RT_List':'RT value in each sample'},inplace=True)
        Filter_RefList_output.rename(columns={'RI_List':'RI value in each sample'},inplace=True)
        Filter_RefList_output.rename(columns={'MS2_Int':'Intensity values of product ions'},inplace=True)
        Filter_RefList_output.rename(columns={'MS2_MZ':'m/z values of product ions'},inplace=True)
        Filter_RefList_output.rename(columns={'SimilarityScore':'Average Similarity Score'},inplace=True)
        Filter_RefList_output.rename(columns={'SC_List':'Similarity Score in each sample'},inplace=True)
        Filter_RefList_output.rename(columns={'max_SampleInt':'Max sample intensity'},inplace=True)
        Filter_RefList_output.rename(columns={'mean_BlankInt':'Average blank intensity'},inplace=True)
        Filter_RefList_output.to_excel(self.filepath_title+'Alignment-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        self.Label_process_sub.setText('Finished')
        self.process_bar.setMaximum(100)
        self.process_bar.setValue(100)
        self.PushButton_Run.setEnabled(True)
        
    
    def Run(self):
        self.PushButton_Run.setEnabled(False)
        self.SampleData = {}
        self.BlankData = {}
        self.SampleClass = {}
        self.MS2_Data = {}
        self.QCData = {}
        for i in self.Data_params.keys():
            if self.Data_params[i]['Type'] == 'Sample':
                self.SampleData[i] = self.Data_params[i]
                self.SampleClass[i] = self.Data_params[i]['Group']
            elif self.Data_params[i]['Type'] == 'MS2':
                self.MS2_Data[i] = self.Data_params[i]
            elif self.Data_params[i]['Type'] == 'Blank':
                self.BlankData[i] = self.Data_params[i]
            elif self.Data_params[i]['Type'] == 'QC':
                self.QCData[i] = self.Data_params[i]
        for i in self.SampleData.keys():
            name_begin = self.SampleData[i]['Path'].rfind('/')
            if name_begin > 0:      
                self.filepath_title = self.SampleData[i]['Path'][0:name_begin]+'/'
            else:
                name_begin = self.SampleData[i]['Path'].rfind('\\')
                self.filepath_title = self.SampleData[i]['Path'][0:name_begin]+'\\'
            break
        self.ClassList = []
        for i in self.SampleClass.keys():
            self.ClassList.append(self.SampleClass[i])
        self.ClassList = list(set(self.ClassList))
        self.Align_path = self.filepath_title+'Alignment-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx'
        ''' --- Heuristic List ---'''
        self.Label_process_sub.setText('Set Heuristic List')
        self.process_bar.setMaximum(100)
        self.process_bar.setValue(0)
        MS_Tor = float(self.LineEdit_MS_Tor.text())/1000000
        with open(self.HeuristicListPath,'rb') as f:
            Alltemp_HPL = pickle.load(f)
        Alltemp_HPL.sort_values('AverageMZ',inplace=True)
        temp_Heuristic_peak_list = []
        for ii in range(len(Alltemp_HPL)):
            self.process_bar.setValue(int(ii/len(Alltemp_HPL))*100)
            QApplication.processEvents()
            if len(temp_Heuristic_peak_list) == 0:
                temp_Heuristic_peak_list = Alltemp_HPL.iloc[[0],:].copy()
                temp_Heuristic_peak_list['PeakNumber'] = 1
            else:
                MZ = Alltemp_HPL.at[ii,'AverageMZ']
                RI_left = Alltemp_HPL.at[ii,'RI_left']
                RI_right = Alltemp_HPL.at[ii,'RI_right']
                RI = Alltemp_HPL.at[ii,'RI']
                MZ_same_list = range(bisect.bisect_left(temp_Heuristic_peak_list['AverageMZ'],MZ*(1-MS_Tor)),bisect.bisect_right(temp_Heuristic_peak_list['AverageMZ'],MZ*(1+MS_Tor)))
                same_list = list(filter(lambda x:((RI_left<=temp_Heuristic_peak_list.at[x,'RI_right'] and RI_right>=temp_Heuristic_peak_list.at[x,'RI_right']) or 
                                         (RI_left<=temp_Heuristic_peak_list.at[x,'RI_left'] and RI_right>=temp_Heuristic_peak_list.at[x,'RI_left'])),MZ_same_list))
                contain_list = list(filter(lambda x:(RI_left>=temp_Heuristic_peak_list.at[x,'RI_left'] and RI_right<=temp_Heuristic_peak_list.at[x,'RI_right']),MZ_same_list))
                if len(same_list) == 0 :
                    if len(contain_list) == 0:
                        iii = Alltemp_HPL.loc[[ii],:].copy()
                        iii['PeakNumber'] = 1
                        temp_Heuristic_peak_list = pd.concat([temp_Heuristic_peak_list,iii],ignore_index=True)
                        temp_Heuristic_peak_list.loc[len(temp_Heuristic_peak_list)-1,'RI_List'] = [RI]
                    elif len(contain_list) == 1:
                        temp_Heuristic_peak_list.loc[contain_list[0],'PeakNumber'] = temp_Heuristic_peak_list.loc[contain_list[0],'PeakNumber']+1
                        temp_Heuristic_peak_list.loc[contain_list[0],'MS_List'].append(MZ)
                        temp_Heuristic_peak_list.loc[contain_list[0],'RI_List'].append(RI)
                        temp_Heuristic_peak_list.loc[contain_list[0],'AverageMZ'] = np.mean(temp_Heuristic_peak_list.loc[contain_list[0],'MS_List'])
                elif len(same_list) == 1 :
                    if RI_left<temp_Heuristic_peak_list.at[same_list[0],'RI_left']:
                        temp_Heuristic_peak_list.loc[same_list[0],'RI_left'] = RI_left
                    if RI_right>temp_Heuristic_peak_list.at[same_list[0],'RI_right']:
                        temp_Heuristic_peak_list.loc[same_list[0],'RI_right'] = RI_right
                    temp_Heuristic_peak_list.loc[same_list[0],'MS_List'].append(MZ)
                    temp_Heuristic_peak_list.loc[same_list[0],'RI_List'].append(RI)
                    temp_Heuristic_peak_list.loc[same_list[0],'AverageMZ'] = np.mean(temp_Heuristic_peak_list.loc[same_list[0],'MS_List'])
                    temp_Heuristic_peak_list.loc[same_list[0],'PeakNumber'] = temp_Heuristic_peak_list.loc[same_list[0],'PeakNumber']+1
                elif len(same_list) == 2 :
                    iii = Alltemp_HPL.loc[[ii],:].copy()
                    iii['PeakNumber'] = temp_Heuristic_peak_list.loc[same_list[0],'PeakNumber'] + temp_Heuristic_peak_list.loc[same_list[1],'PeakNumber'] + 1
                    iii.loc[ii,'RI_right'] = max([RI_right,temp_Heuristic_peak_list.loc[same_list[0],'RI_right'],temp_Heuristic_peak_list.loc[same_list[1],'RI_right']])
                    iii.loc[ii,'RI_left'] = min([RI_left,temp_Heuristic_peak_list.loc[same_list[0],'RI_left'],temp_Heuristic_peak_list.loc[same_list[1],'RI_left']])
                    iii['MS_List'] = iii['MS_List'].astype('object')
                    iii.at[ii,'MS_List'] = list(temp_Heuristic_peak_list.loc[same_list[0],'MS_List'] + temp_Heuristic_peak_list.loc[same_list[1],'MS_List'] + iii.loc[ii,'MS_List'])
                    iii['RI_List'] = iii['RI_List'].astype('object')
                    iii.at[ii,'RI_List'] = list(temp_Heuristic_peak_list.loc[same_list[0],'RI_List'] + temp_Heuristic_peak_list.loc[same_list[1],'RI_List'] + [iii.loc[ii,'RI']])
                    iii.loc[ii,'AverageMZ'] = float(np.mean(iii.loc[ii,'MS_List']))
                    temp_Heuristic_peak_list.drop(same_list,inplace=True)
                    temp_Heuristic_peak_list = pd.concat([temp_Heuristic_peak_list,iii],ignore_index=True)
        self.HeuristicList =  temp_Heuristic_peak_list.copy()
        ''' --- Process Sample Data ---'''
        self.Label_process_sub.setText('Load Data')
        self.process_bar.setMaximum(len(self.SampleData.keys()))
        self.process_bar.setValue(0)
        
        self.params = {
            'Path':self.SampleData,
            'MS1_Tor':float(self.LineEdit_MS_Tor.text())/1000000,
            'RI_Tor':float(self.LineEdit_RI_Tor.text())/100,
            'RT_Tor':float(self.LineEdit_RT_Tor.text())*60,
            'min_Int':int(self.LineEdit_Int_min.text()),
            'Points':int(self.LineEdit_Point.text()),
            'RIISPath':self.RIISPath,
            'HeuristicList':self.HeuristicList.copy(),
            'SN':int(self.LineEdit_SN.text()),
        }
        self.worker_sample = PeakDetectWorker(self.SampleData, self.params)
        self.worker_sample.progress.connect(self.update_progress)
        self.worker_sample.finished.connect(self.on_PeakDetectWorker_done)
        self.worker_sample.start()
        
        
class PeakDetectWorker(QThread):
    progress = pyqtSignal(int)       # 当前进度
    finished = pyqtSignal(dict)      # 所有结果完成

    def __init__(self, Data, params, parent=None):
        super().__init__(parent)
        self.Data = Data
        self.params = params

    def run(self):
        total_tasks = len(self.Data)
        result_dict = {}
        # 使用 ProcessPoolExecutor 提交任务
        with ProcessPoolExecutor(max_workers=max(1, os.cpu_count()-1)) as executor:
            future_to_key = {}
            for key, info in self.Data.items():
                future = executor.submit(
                    pool_HP,
                    info['Path'],
                    self.params['MS1_Tor'],
                    self.params['RI_Tor'],
                    self.params['RT_Tor'],
                    self.params['min_Int'],
                    self.params['Points'],
                    self.params['RIISPath'],
                    self.params['HeuristicList'],
                    self.params['SN']
                )
                future_to_key[future] = key
            tasks_done = 0
            # as_completed 会按完成顺序迭代返回结果
            for future in as_completed(future_to_key):
                key = future_to_key[future]
                try:
                    res = future.result()  # 返回 tuple of DataFrames
                    # 用 copy() 避免 Manager/dict pickle 问题
                    result_dict[key] = (res[0].copy(), res[1].copy())
                except Exception as e:
                    print(f"Error in {key}: {e}")
                tasks_done += 1
                percent = int(tasks_done / total_tasks * 100)
                self.progress.emit(percent)
        # 所有任务完成，发射 finished 信号
        self.finished.emit(result_dict)

class MS2Worker(QThread):
    progress = pyqtSignal(int)       # 当前进度
    finished = pyqtSignal(dict)      # 所有结果完成
    def __init__(self, Data, params, parent=None):
        super().__init__(parent)
        self.Data = Data
        self.params = params
    def run(self):
        total_tasks = len(self.Data)
        result_dict = {}
        with ProcessPoolExecutor(max_workers=max(1, os.cpu_count()-1)) as executor:
            future_to_key = {}
            for key, info in self.Data.items():
                future = executor.submit(
                    pool_load_MS2_Data,
                    info['Path'],
                    self.params['MS1_Tor'],
                    self.params['RI_Tor'],
                    self.params['RT_Tor'],
                    self.params['min_Int'],
                    self.params['Points'],
                    self.params['RIISPath']
                )
                future_to_key[future] = key
            tasks_done = 0
            # as_completed 会按完成顺序迭代返回结果
            for future in as_completed(future_to_key):
                key = future_to_key[future]
                try:
                    res = future.result()  # 返回 tuple of DataFrames
                    # 用 copy() 避免 Manager/dict pickle 问题
                    result_dict[key] = res.copy()
                except Exception as e:
                    print(f"Error in {key}: {e}")
                tasks_done += 1
                percent = int(tasks_done / total_tasks * 100)
                self.progress.emit(percent)
        # 所有任务完成，发射 finished 信号
        self.finished.emit(result_dict)



class ProgressBar(QDialog):
    def __init__(self,parent=None):
        super(ProgressBar,self).__init__(parent)
        self.resize(500,32)
        self.progressBar = QProgressBar(self)
        self.progressBar.setMaximum(500)
        self.progressBar.setMinimum(0)
        self.progressBar.setValue(0)
        self.centerWindow()
        self.show()
        
    def setValue(self,task_Number,total_task_Number,value):
        if total_task_Number == 1:
            self.setWindowTitle('Processing')
        else:
            self.setWindowTitle('Processing '+str(task_Number)+'/'+str(total_task_Number))
        self.progressBar.setValue(value)
        
    def centerWindow(self):
        screen = QDesktopWidget().screenGeometry()
        size  = self.geometry()
        LeftValue  = int((screen.width()-size.width())/2)
        TopValue = int((screen.height()-size.height())/2)
        self.move(LeftValue,TopValue)

def pool_HP(Path,MS1_Tor,RI_Tor,RT_Tor,min_Int,Points,RIISPath,HeuristicList,SN):
    temp_EMZDP = EazyMZDataProcess(Path)
    temp_EMZDP.set_param('MS1_Tor',MS1_Tor)
    temp_EMZDP.set_param('RI_Tor',RI_Tor)
    temp_EMZDP.set_param('RT_Tor',RT_Tor)
    temp_EMZDP.set_param('min_Int',min_Int)
    temp_EMZDP.set_param('Points',Points)
    temp_EMZDP.set_RIIS(RIISPath)
    File_path = glob.glob(Path.replace('.mzML','.pkl'))
    if len(File_path) > 0:
        temp_EMZDP.load_FPD(File_path[0])
    else:
        temp_EMZDP.Heuristic_peak_list = HeuristicList
        temp_EMZDP.Final_Peak_Detect = temp_EMZDP.Heuristic_PeakDetect()
        temp_EMZDP.Calculate_RI()
        temp_EMZDP.Calculate_SN(drop=True,Threshold=SN)
        temp_EMZDP.save_FPD(Path.replace('.mzML','.pkl'))
    return temp_EMZDP.Final_Peak_Detect.copy(),temp_EMZDP.RIIS.copy()

def pool_load_MS2_Data(Path,MS1_Tor,RI_Tor,RT_Tor,min_Int,Points,RIISPath):
    temp_EMZDP = EazyMZDataProcess(Path)
    temp_EMZDP.set_param('MS1_Tor',MS1_Tor)
    temp_EMZDP.set_param('RI_Tor',RI_Tor)
    temp_EMZDP.set_param('RT_Tor',RT_Tor)
    temp_EMZDP.set_param('min_Int',min_Int)
    temp_EMZDP.set_param('Points',Points)
    temp_EMZDP.set_RIIS(RIISPath)
    temp_EMZDP.Calculate_MS2_RI()
    return temp_EMZDP.MS2_Data.copy()

class EazyMZDataProcess(object):   
    def __init__(self,DataPath):
        if str(DataPath.rfind('/')) != -1 :
            self.DataName = DataPath[str(DataPath).rfind('/')+1:len(str(DataPath))-5]
        else:
            self.DataName = DataPath[str(DataPath).rfind('\\')+1:len(str(DataPath))-5]
        self.OriginData = pyopenms.MSExperiment()
        self.file_path = DataPath
        ''' 储存文件pyopenms.MzMLFile().store("filtered.mzML", exp) '''
        if self.file_path.endswith('mzML'):
            pyopenms.MzMLFile().load(self.file_path,self.OriginData)
        elif self.file_path.endswith('mzXML'):
            pyopenms.MzXMLFile().load(self.file_path,self.OriginData)
        self.OriginData.sortSpectra(True)
        self.__param = {'MS1_Tor':0.000010,'RT_Tor':6,'min_Int':10000,'min_RT':0,'min_RT_width':6,'max_Noise':2000,
                        'Deconvolution':False,'FeatureDetectPlot':3,'MergeRule':'Intersection',
                        'UpDown_gap':10,'saveAutoList':False,'smooth':5,'Points':25,'RI_Tor':0.02}
        self.Origin_RT_List = np.array([])
        self.Origin_MZ_List = []
        self.Origin_Int_List = []
        self.MS2_Pre = []
        self.MS2_RT_List = []
        self.MS2_MZ_List = []
        self.MS2_Int_List = []
        self.MS2_RelInt=[]
        for i in self.OriginData:
            if i.getMSLevel()==1:
                MZ_temp, Int_temp = i.get_peaks()
                self.Origin_MZ_List.append(np.around(MZ_temp,5))
                self.Origin_Int_List.append(np.around(Int_temp,0))
                self.Origin_RT_List = np.append(self.Origin_RT_List,i.getRT())
            if i.getMSLevel()==2:
                Pre_temp = i.getPrecursors()[0].getMZ() 
                MZ_temp,Int_temp=i.get_peaks()
                temp_range = range(1,bisect.bisect(MZ_temp,Pre_temp*(1+self.get_param('MS1_Tor'))))
                temp_range = list(filter(lambda x:abs(MZ_temp[x-1]-MZ_temp[x])/MZ_temp[x] > self.get_param('MS1_Tor') or Int_temp[x-1] < Int_temp[x],temp_range))
                if len(temp_range)>0:
                    if abs(MZ_temp[0]-MZ_temp[temp_range[0]])/MZ_temp[temp_range[0]] > self.get_param('MS1_Tor') or Int_temp[0] > Int_temp[temp_range[0]]:
                        temp_range = [0]+temp_range
                    MZ_temp = MZ_temp[temp_range]
                    Int_temp = Int_temp[temp_range]
                    self.MS2_Pre.append(Pre_temp)
                    self.MS2_RT_List.append(i.getRT())
                    self.MS2_MZ_List.append(MZ_temp)
                    self.MS2_Int_List.append(Int_temp)
                    self.MS2_RelInt.append(max(Int_temp))
                '''
                else:
                    MZ_temp = []
                    Int_temp = []
                    self.MS2_Pre.append(Pre_temp)
                    self.MS2_RT_List.append(i.getRT())
                    self.MS2_MZ_List.append(MZ_temp)
                    self.MS2_Int_List.append(Int_temp)
                    self.MS2_RelInt.append(0)
                '''
        self.Origin_RT_List = np.around(self.Origin_RT_List,3)
        self.MS2_RT_List = np.around(self.MS2_RT_List,3)
        self.MS1_Data = {'Scan_Time':self.Origin_RT_List,'MZ_List':self.Origin_MZ_List,'Int_List':self.Origin_Int_List}
        self.MS1_Data = pd.DataFrame(self.MS1_Data)
        self.MS2_Data = {'Scan_Time':self.MS2_RT_List,'Pre_MZ':self.MS2_Pre,'MZ_List':self.MS2_MZ_List,'Int_List':self.MS2_Int_List,'Rel_Int':self.MS2_RelInt}
        self.MS2_Data = pd.DataFrame(self.MS2_Data)
        self.__param['Flow_RT'] = int(self.__param['min_RT_width']/(4*sum(np.diff(self.Origin_RT_List))/len(np.diff(self.Origin_RT_List))))+1
        if self.OriginData[0].getInstrumentSettings().getPolarity() == 1:
            self.__param['Polarity'] = 'Positive'
        elif self.OriginData[0].getInstrumentSettings().getPolarity() == 2:
            self.__param['Polarity'] = 'Negative'
        else:
            self.__param['Polarity'] = 'Not give'
            print('Uncertain Polarity')
        self.OriginData = []
    def add_0(x):
        x.append(0)
        return x
    def add_Blank(self,DataPath):
        self.BlankData = pyopenms.MSExperiment()
        if DataPath.endswith('mzML'):
            pyopenms.MzMLFile().load(DataPath,self.BlankData)
        elif DataPath.endswith('mzXML'):
            pyopenms.MzXMLFile().load(DataPath,self.BlankData)
        self.BlankData.sortSpectra(True)
        self.Blank_RT_List = np.array([])
        self.Blank_MZ_List = []
        self.Blank_Int_List = []
        for i in self.BlankData:
            if i.getMSLevel()==1:
                MZ_temp, Int_temp = i.get_peaks()
                self.Blank_MZ_List.append(np.around(MZ_temp,5))
                self.Blank_Int_List.append(np.around(Int_temp,0))
                self.Blank_RT_List = np.append(self.Blank_RT_List,i.getRT())
        self.Blank_RT_List = np.around(self.Blank_RT_List,3)
    def add_RT(x, RT):
        x.append(RT)
        return x
    def add_param(self,name,value):
        self.__param[name]=value
        print('set',name,' = ',value)
    def add_Heuristic_peak_list(self,path):
        Heuristic_peak_list = pd.read_excel(path)
        Heuristic_peak_list['RI_left'] = Heuristic_peak_list['RI'].apply(lambda x:x*(1-self.get_param('RI_Tor')))
        Heuristic_peak_list['RI_right'] = Heuristic_peak_list['RI'].apply(lambda x:x*(1+self.get_param('RI_Tor')))
        Heuristic_peak_list['MZ_List'] = Heuristic_peak_list['AverageMZ'].apply(lambda x:[x])
        Heuristic_peak_list['RI_List'] = Heuristic_peak_list['RI'].apply(lambda x:[x])
        Alltemp_HPL=pd.DataFrame(columns=['AverageMZ','RT','MZ_List','RT_List'])
        Alltemp_HPL['AverageMZ'] = Alltemp_HPL['AverageMZ'].map(lambda x:'%.4f'%x)
        sub_HPL_dict = {}
        for i in list(set(Heuristic_peak_list['Gradient Code'])):
            sub_Heuristic_peak_list = Heuristic_peak_list[Heuristic_peak_list['Gradient Code']==i].copy()
            sub_Heuristic_peak_list.reset_index(drop=True,inplace=True)
            sub_HPL_dict[i] = sub_Heuristic_peak_list.copy()
        for i in sub_HPL_dict.keys():
            sub_HPL_dict[i].sort_values('AverageMZ',ascending=(False),inplace=True)
            sub_HPL_dict[i].reset_index(drop=True,inplace=True)
            Alltemp_HPL.sort_values('AverageMZ',ascending=(False),inplace=True)
            Alltemp_HPL.reset_index(drop=True,inplace=True)
            temp_Data = sub_HPL_dict[i].copy()
            init_RefMZ = pd.DataFrame(columns=['AverageMZ','MZList','RowIndex'])
            init_SampleMZ = pd.DataFrame(columns=['AverageMZ','MZList','RowIndex'])
            Sample_Time = temp_Data.loc[:,'RI']
            for ii in range(len(Alltemp_HPL)):
                if len(init_RefMZ)==0:
                    init_RefMZ.loc[len(init_RefMZ)] = [Alltemp_HPL.at[ii,'AverageMZ'],[Alltemp_HPL.at[ii,'AverageMZ']],[ii]]
                else:
                    if abs(init_RefMZ.at[len(init_RefMZ)-1,'AverageMZ']-Alltemp_HPL.at[ii,'AverageMZ'])/Alltemp_HPL.at[ii,'AverageMZ']<self.__param['MS1_Tor']:
                        init_RefMZ.at[len(init_RefMZ)-1,'MZList'].append(Alltemp_HPL.at[ii,'AverageMZ'])
                        init_RefMZ.at[len(init_RefMZ)-1,'RowIndex'].append(ii)
                        init_RefMZ.at[len(init_RefMZ)-1,'AverageMZ']=np.mean(init_RefMZ.at[len(init_RefMZ)-1,'MZList'])
                    else:
                        init_RefMZ.loc[len(init_RefMZ)] = [Alltemp_HPL.at[ii,'AverageMZ'],[Alltemp_HPL.at[ii,'AverageMZ']],[ii]]
            for ii in range(len(temp_Data)):
                if len(init_SampleMZ)==0:
                    init_SampleMZ.loc[len(init_SampleMZ)] = [temp_Data.at[ii,'AverageMZ'],[temp_Data.at[ii,'AverageMZ']],[ii]]
                else:
                    if abs(init_SampleMZ.at[len(init_SampleMZ)-1,'AverageMZ']-temp_Data.at[ii,'AverageMZ'])/temp_Data.at[ii,'AverageMZ']<self.__param['MS1_Tor']:
                        init_SampleMZ.at[len(init_SampleMZ)-1,'MZList'].append(temp_Data.at[ii,'AverageMZ'])
                        init_SampleMZ.at[len(init_SampleMZ)-1,'RowIndex'].append(ii)
                        init_SampleMZ.at[len(init_SampleMZ)-1,'AverageMZ']=np.mean(init_SampleMZ.at[len(init_SampleMZ)-1,'MZList'])
                    else:
                        init_SampleMZ.loc[len(init_SampleMZ)] = [temp_Data.at[ii,'AverageMZ'],[temp_Data.at[ii,'AverageMZ']],[ii]]
            add_MZ = []
            add_RT = []
            add_MS_List = []
            add_RT_List = []  
            for ii in range(len(init_SampleMZ)):
                Sample_MZ = init_SampleMZ.at[ii,'AverageMZ']
                match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'AverageMZ']-Sample_MZ)/Sample_MZ<self.__param['MS1_Tor'],range(len(init_RefMZ))))
                if len(match_Index)>1:
                    Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                    for i_mI in range(1,len(match_Index)):
                        Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                elif len(match_Index)==1:
                    Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                else:
                    for i_add in init_SampleMZ.at[ii,'RowIndex']:
                        add_MZ.append(temp_Data.at[i_add,'AverageMZ'])
                        
                        add_RT.append(Sample_Time[i_add])
                        add_MS_List.append([temp_Data.at[i_add,'AverageMZ']])
                        add_RT_List.append([Sample_Time[i_add]])
                    continue
                Score_matrix= np.zeros([len(Ref_Index),len(init_SampleMZ.at[ii,'RowIndex'])])
                for i_Ref in range(len(Ref_Index)):
                    for i_Sample in range(len(init_SampleMZ.at[ii,'RowIndex'])):
                        i_Ref_MZ = Alltemp_HPL.at[Ref_Index[i_Ref],'AverageMZ']
                        i_Sample_MZ = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sample],'AverageMZ']
                        i_Ref_Time = Alltemp_HPL.at[Ref_Index[i_Ref],'RT']
                        i_Sample_Time = Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sample]]
                        if abs(i_Ref_MZ-i_Sample_MZ)/i_Ref_MZ<self.get_param('MS1_Tor') and abs(i_Ref_Time-i_Sample_Time)/i_Ref_Time<self.get_param('RI_Tor'):
                            Score_matrix[i_Ref,i_Sample] = 0.5*np.exp(-0.5*((i_Sample_Time-i_Ref_Time)/(i_Ref_Time*self.get_param('RI_Tor')))**2)+(1-0.5)*np.exp(-0.5*((i_Sample_MZ-i_Ref_MZ)/(i_Sample_MZ*self.__param['MS1_Tor']))**2)
                Sm_row,Sm_col = linear_sum_assignment(Score_matrix,True)
                for i_Sm_row,i_Sm_col in zip(Sm_row,Sm_col):
                    if Score_matrix[i_Sm_row,i_Sm_col] != 0:
                        Alltemp_HPL.at[Ref_Index[i_Sm_row],'MZ_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                        Alltemp_HPL.at[Ref_Index[i_Sm_row],'RT_List'].append(Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                        Alltemp_HPL.at[Ref_Index[i_Sm_row],'RT'] = np.mean(Alltemp_HPL.at[Ref_Index[i_Sm_row],'RT_List'])
                        Alltemp_HPL.at[Ref_Index[i_Sm_row],'AverageMZ'] = np.mean(Alltemp_HPL.at[Ref_Index[i_Sm_row],'MZ_List'])
                    else:
                        add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                        
                        add_RT.append(Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                        add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ']])
                        add_RT_List.append([Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]]])
                miss_col = list(filter(lambda x:x not in Sm_col,range(len(init_SampleMZ.at[ii,'RowIndex']))))
                for i_miss in miss_col:
                    add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ'])
                    
                    add_RT.append(Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_miss]])
                    add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ']])
                    add_RT_List.append([Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_miss]]])
            add_RefList = pd.DataFrame({'AverageMZ': add_MZ,'RT': add_RT,'MZ_List':add_MS_List,'RT_List':add_RT_List})
            Alltemp_HPL = pd.concat([Alltemp_HPL, add_RefList])
            Alltemp_HPL.reset_index(drop=True, inplace=True)
            Alltemp_HPL = Alltemp_HPL.fillna(0)
        Alltemp_HPL['RT'] = Alltemp_HPL['RT_List'].apply(lambda x:sum(x)/len(x))
        Alltemp_HPL.rename(columns={'RT':'RI'},inplace=True)
        Alltemp_HPL.rename(columns={'RT_List':'RI_List'},inplace=True)
        Alltemp_HPL['RI_left'] = Alltemp_HPL['RI'].apply(lambda x:x*(1-self.get_param('RI_Tor')))
        Alltemp_HPL['RI_right'] = Alltemp_HPL['RI'].apply(lambda x:x*(1+self.get_param('RI_Tor')))
        Alltemp_HPL.sort_values('AverageMZ',ascending=(True),inplace=True)
        Alltemp_HPL.reset_index(drop=True,inplace=True)
        self.Heuristic_normal_list = Alltemp_HPL.copy()
        temp_Heuristic_peak_list = []
        for ii in range(len(Alltemp_HPL)):
            if len(temp_Heuristic_peak_list) == 0:
                temp_Heuristic_peak_list = Alltemp_HPL.iloc[[0],:].copy()
                temp_Heuristic_peak_list['PeakNumber'] = 1
            else:
                MZ = Alltemp_HPL.at[ii,'AverageMZ']
                RI_left = Alltemp_HPL.at[ii,'RI_left']
                RI_right = Alltemp_HPL.at[ii,'RI_right']
                RI = Alltemp_HPL.at[ii,'RI']
                same_list = list(filter(lambda x:EazyMZDataProcess.ppm_compare(MZ,temp_Heuristic_peak_list.at[x,'AverageMZ'])<self.__param['MS1_Tor'] and 
                                        ((RI_left<=temp_Heuristic_peak_list.at[x,'RI_right'] and RI_right>=temp_Heuristic_peak_list.at[x,'RI_right']) or 
                                         (RI_left<=temp_Heuristic_peak_list.at[x,'RI_left'] and RI_right>=temp_Heuristic_peak_list.at[x,'RI_left'])),range(len(temp_Heuristic_peak_list))))
                contain_list = list(filter(lambda x:EazyMZDataProcess.ppm_compare(MZ,temp_Heuristic_peak_list.at[x,'AverageMZ'])<self.__param['MS1_Tor'] and (RI_left>=temp_Heuristic_peak_list.at[x,'RI_left'] and RI_right<=temp_Heuristic_peak_list.at[x,'RI_right']),range(len(temp_Heuristic_peak_list))))
                if len(same_list) == 0 :
                    if len(contain_list) == 0:
                        iii = Alltemp_HPL.iloc[[ii],:].copy()
                        iii['PeakNumber'] = 1
                        temp_Heuristic_peak_list = pd.concat([temp_Heuristic_peak_list,iii])
                        temp_Heuristic_peak_list.reset_index(drop=True,inplace=True)
                        temp_Heuristic_peak_list.loc[len(temp_Heuristic_peak_list)-1,'RI_List'] = [RI]
                    elif len(contain_list) == 1:
                        temp_Heuristic_peak_list.loc[contain_list[0],'PeakNumber'] = temp_Heuristic_peak_list.loc[contain_list[0],'PeakNumber']+1
                        temp_Heuristic_peak_list.loc[contain_list[0],'MZ_List'].append(MZ)
                        temp_Heuristic_peak_list.loc[contain_list[0],'RI_List'].append(RI)
                        temp_Heuristic_peak_list.loc[contain_list[0],'AverageMZ'] = np.mean(temp_Heuristic_peak_list.loc[contain_list[0],'MZ_List'])
                elif len(same_list) == 1 :
                    if RI_left<temp_Heuristic_peak_list.at[same_list[0],'RI_left']:
                        temp_Heuristic_peak_list.loc[same_list[0],'RI_left'] = RI_left
                    if RI_right>temp_Heuristic_peak_list.at[same_list[0],'RI_right']:
                        temp_Heuristic_peak_list.loc[same_list[0],'RI_right'] = RI_right
                    temp_Heuristic_peak_list.loc[same_list[0],'MZ_List'].append(MZ)
                    temp_Heuristic_peak_list.loc[same_list[0],'RI_List'].append(RI)
                    temp_Heuristic_peak_list.loc[same_list[0],'AverageMZ'] = np.mean(temp_Heuristic_peak_list.loc[same_list[0],'MZ_List'])
                    temp_Heuristic_peak_list.loc[same_list[0],'PeakNumber'] = temp_Heuristic_peak_list.loc[same_list[0],'PeakNumber']+1
                elif len(same_list) == 2 :
                    iii = Alltemp_HPL.iloc[[ii],:].copy()
                    iii['PeakNumber'] = temp_Heuristic_peak_list.loc[same_list[0],'PeakNumber'] + temp_Heuristic_peak_list.loc[same_list[1],'PeakNumber'] + 1
                    iii.loc[ii,'RI_right'] = max([RI_right,temp_Heuristic_peak_list.loc[same_list[0],'RI_right'],temp_Heuristic_peak_list.loc[same_list[1],'RI_right']])
                    iii.loc[ii,'RI_left'] = min([RI_left,temp_Heuristic_peak_list.loc[same_list[0],'RI_left'],temp_Heuristic_peak_list.loc[same_list[1],'RI_left']])
                    iii['MZ_List'] = iii['MZ_List'].astype('object')
                    iii.at[ii,'MZ_List'] = list(temp_Heuristic_peak_list.loc[same_list[0],'MZ_List'] + temp_Heuristic_peak_list.loc[same_list[1],'MZ_List'] + iii.loc[ii,'MZ_List'])
                    iii['RI_List'] = iii['RI_List'].astype('object')
                    iii.at[ii,'RI_List'] = list(temp_Heuristic_peak_list.loc[same_list[0],'RI_List'] + temp_Heuristic_peak_list.loc[same_list[1],'RI_List'] + [iii.loc[ii,'RI']])
                    iii.loc[ii,'AverageMZ'] = float(np.mean(iii.loc[ii,'MZ_List']))
                    temp_Heuristic_peak_list.drop(same_list,inplace=True)
                    temp_Heuristic_peak_list = pd.concat([temp_Heuristic_peak_list,iii])
                    temp_Heuristic_peak_list.reset_index(drop=True,inplace=True)
        self.Heuristic_peak_list = temp_Heuristic_peak_list.copy()
    def ThreeZero(x):
        if len(x) >= 3:
            if x[-1] == 0 and x[-2] == 0 and x[-3] == 0:
                return False
            else:
                return True
        else:
            return True
    def ClosestPosition(TargetNumber, List):
        if TargetNumber >= List[-1]:
            #print("Close error +")
            return len(List)-1
        elif TargetNumber <= List[0]:
            #print("Close error -")
            return 0
        position = bisect.bisect_left(List, TargetNumber)
        before = List[position-1]
        after = List[position]
        if after - TargetNumber < TargetNumber - before:
            return position
        else:
            return position-1
    def ppm_compare(a, b):
        ppm_result = abs(a-b)/b
        return ppm_result
    '''
    def Draw_TIC(self):
        TIC_Int_List = []
        for i in range(len(self.Origin_Int_List)):
            TIC_Int_List.append(sum(self.Origin_Int_List[i]))
        plt.plot(self.Origin_RT_List,TIC_Int_List)
        plt.show()
        TIC_Table = {'ScanTime':self.Origin_RT_List,'Intensity':TIC_Int_List}
        TIC_Table = pd.DataFrame(TIC_Table)
        #TIC_Table.to_excel(OutputPath+'TIC.xlsx')
        return TIC_Table
    '''
    def CosineSimilarity(ExestList, NewList, LOG=False):
        if  len(ExestList) == 0 :
            return 0
        else:
            if len(ExestList) != len(NewList):
                print("CosineSimilarity error", len(ExestList), len(NewList))
                return 0
            if LOG == True:
                temp_ExestList = []
                for i in range(len(ExestList)):
                    if ExestList[i] != 0:
                        temp_ExestList.append(math.log(ExestList[i]))
                    else:
                        temp_ExestList.append(0)
                ExestList = temp_ExestList
                temp_NewList = []
                for ii in range(len(NewList)):
                    if NewList[ii] != 0:
                        temp_NewList.append(math.log(NewList[ii]))
                    else:
                        temp_NewList.append(0)
                NewList = temp_NewList
            ExestList = np.array(ExestList)
            NewList = np.array(NewList)
            CosineSimilarityValue = ExestList.dot(
                NewList)/(np.linalg.norm(ExestList) * np.linalg.norm(NewList))
            return CosineSimilarityValue
    def Calculate_SB(self,drop=False,limit=5):
        self.Final_Peak_Detect['BLANK'] = 0
        for i in range(len(self.Final_Peak_Detect)):
            RTL = self.Final_Peak_Detect.at[i,'RTList'][0]
            RTR = self.Final_Peak_Detect.at[i,'RTList'][-1]
            MZ = self.Final_Peak_Detect.at[i,'AverageMZ']
            [RT_List,Int_List] = self.ExtractBlankPoint(MZ,RTR,RTL,smooth_index=5)
            self.Final_Peak_Detect.at[i,'BLANK'] = max(1,max(Int_List))
        self.Final_Peak_Detect['S/B'] = self.Final_Peak_Detect['Int']/self.Final_Peak_Detect['BLANK']
        if drop == True:
            self.Final_Peak_Detect = self.Final_Peak_Detect[self.Final_Peak_Detect['S/B']>limit].copy()
            self.Final_Peak_Detect.reset_index(drop=True,inplace=True)
    def ExtractDataPoint(self,MZ, RTR, RTL,s_min=False, plt_for_test=False,smooth_index=5):
        smooth_index =(smooth_index-1)//2
        if smooth_index<0:
            smooth_index=0
        if s_min == True:
            RTR = RTR * 60          # min -> s
            RTL = RTL * 60
        RTR_place = EazyMZDataProcess.ClosestPosition(RTR,self.Origin_RT_List)
        RTL_place = EazyMZDataProcess.ClosestPosition(RTL,self.Origin_RT_List)
        RT_List = self.Origin_RT_List[RTL_place:RTR_place]
        Int_List = []
        #MZ_List=[]
        for i in range(RTL_place, RTR_place):
            if len(self.Origin_MZ_List[i])>0:
                if len(self.Origin_MZ_List[i])>=4:
                    #MZ_place = EazyMZDataProcess.ClosestPosition(MZ,self.Origin_MZ_List[i])
                    MZ_place = bisect.bisect_right(self.Origin_MZ_List[i],MZ)
                    if MZ_place > len(self.Origin_MZ_List[i])-3 and MZ_place>=3:
                        MZ_place = len(self.Origin_MZ_List[i])-3
                    elif MZ_place < 2:
                        MZ_place = 2
                    temp_int = []
                    #temp_MZ = []
                    for ii in range(MZ_place-2, MZ_place+2):
                        if MZ*(1-self.__param['MS1_Tor']) < self.Origin_MZ_List[i][ii] < MZ*(1+self.__param['MS1_Tor']):
                            temp_int.append(self.Origin_Int_List[i][ii])
                            #temp_MZ.append(Origin_MZ_List[i][ii])
                        else:
                            temp_int.append(0)
                    if not temp_int:
                        Int_List.append(0)
                        #MZ_List.append(0)
                    else:
                        #temp_int = np.array(temp_int)
                        Int_List.append(max(temp_int))
                        #MZ_List.append()
                else:
                    Int_List.append(0)
            else:
                Int_List.append(0)
        # smooth
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        if smooth_index != 0:
            Int_List = list(map(lambda x:GaussSmooth(Int_List[x-smooth_index:x+smooth_index+1]) if x in range(smooth_index,len(Int_List)-smooth_index) else Int_List[x],range(len(Int_List))))
        '''
        if plt_for_test == True:
            plt.figure()
            plt.title(str(MZ))
            plt.plot(RT_List,Int_List)
            plt.figure()
        '''
        return RT_List, Int_List
    def ExtractBlankPoint(self,MZ, RTR, RTL,s_min=False, plt_for_test=False,smooth_index=5):
        smooth_index =(smooth_index-1)//2
        if smooth_index<0:
            smooth_index=0
        if s_min == True:
            RTR = RTR * 60          # min -> s
            RTL = RTL * 60
        RTR_place = EazyMZDataProcess.ClosestPosition(RTR,self.Blank_RT_List)
        RTL_place = EazyMZDataProcess.ClosestPosition(RTL,self.Blank_RT_List)
        RT_List = self.Blank_RT_List[RTL_place:RTR_place]
        Int_List = []
        #MZ_List=[]
        for i in range(RTL_place, RTR_place):
            if len(self.Blank_MZ_List[i])>0:
                if len(self.Blank_MZ_List[i])>=4:
                    MZ_place = EazyMZDataProcess.ClosestPosition(MZ,self.Blank_MZ_List[i])
                    if MZ_place > len(self.Blank_MZ_List[i])-3 and MZ_place>=3:
                        MZ_place = len(self.Blank_MZ_List[i])-3
                    elif MZ_place < 2:
                        MZ_place = 2
                    temp_int = []
                    #temp_MZ = []
                    for ii in range(MZ_place-2, MZ_place+2):
                        if MZ*(1-self.__param['MS1_Tor']) < self.Blank_MZ_List[i][ii] < MZ*(1+self.__param['MS1_Tor']):
                            temp_int.append(self.Blank_Int_List[i][ii])
                            #temp_MZ.append(Origin_MZ_List[i][ii])
                        else:
                            temp_int.append(0)
                    if not temp_int:
                        Int_List.append(0)
                    else:
                        Int_List.append(max(temp_int))
                else:
                    Int_List.append(0)
            else:
                Int_List.append(0)
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        if smooth_index != 0:
            Int_List = list(map(lambda x:GaussSmooth(Int_List[x-smooth_index:x+smooth_index+1]) if x in range(smooth_index,len(Int_List)-smooth_index) else Int_List[x],range(len(Int_List))))
        '''
        if plt_for_test == True:
            plt.figure()
            plt.title(str(MZ))
            plt.plot(RT_List,Int_List)
            plt.figure()
        '''
        return RT_List, Int_List
    
    def PeakDetecter_for_Heuristic(self,Auto_RT_List,Auto_Int_List,MZ,WhetherAppend=True):
        # 单个色谱峰识别模块，由下方另一函数调用
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        #Final_RTL = Auto_RT_List[0]
        #Final_RTR = Auto_RT_List[-1]
        #Auto_MZ_List = UnprocessList['MZ'][ListIndex]
        if Auto_RT_List[-1]-Auto_RT_List[0] > self.get_param('min_RT_width')*0.8:
            if np.max(Auto_Int_List) > self.get_param('min_Int'):
                Diff_List = np.diff(Auto_Int_List)
                if Diff_List.max() > self.get_param('min_Int')*0.1 and Diff_List.min() < self.get_param('min_Int')*(-0.1):
                    if WhetherAppend == True:
                        RTL = Auto_RT_List[0]-self.get_param('RT_Tor')
                        RTR = Auto_RT_List[-1]+self.get_param('RT_Tor')
                    else:
                        RTL = Auto_RT_List[0]
                        RTR = Auto_RT_List[-1]
                    [Auto_RT_List, Origin_Int_List] = EazyMZDataProcess.ExtractDataPoint(self,MZ, RTR, RTL,s_min=False, plt_for_test=False, smooth_index=0)
                    if self.get_param('smooth')>0:
                        smooth_index =(self.get_param('smooth')-1)//2
                        Auto_Int_List = list(map(lambda x:GaussSmooth(Origin_Int_List[x-smooth_index:x+smooth_index]) if smooth_index<=x<=len(Origin_Int_List)-smooth_index else Origin_Int_List[x],range(len(Origin_Int_List))))
                    else:
                        Auto_Int_List = Origin_Int_List
                    Auto_Int_List = np.array(Auto_Int_List)
                    Auto_Int_List[2:-2] = np.convolve(Auto_Int_List, np.ones(5)/5, mode='valid')
                    Auto_Int_List = list(Auto_Int_List)
                    Diff_List = np.diff(Auto_Int_List[1:len(Auto_Int_List)-2])
                    Diff_List = np.array(list(map(lambda x:GaussSmooth(Diff_List[x-2:x+2+1]) if x in range(2,len(Diff_List)-2) else Diff_List[x],range(len(Diff_List)))))
                    FD_List = np.array(list(map(lambda x: EazyMZDataProcess.TFFD(x, Auto_Int_List, Auto_RT_List), range(2, len(Auto_Int_List)-2))))
                    FD_List = np.array(list(map(lambda x:GaussSmooth(FD_List[x-2:x+2+1]) if x in range(2,len(FD_List)-2) else FD_List[x],range(len(FD_List)))))
                    SD_List = np.array(list(map(lambda x: EazyMZDataProcess.TFSD(x, Auto_Int_List, Auto_RT_List), range(2, len(Auto_Int_List)-2))))
                    SD_List = np.array(list(map(lambda x:GaussSmooth(SD_List[x-2:x+2+1]) if x in range(2,len(SD_List)-2) else SD_List[x],range(len(SD_List)))))
                    ABS_FD_List = abs(np.array(FD_List))
                    FD_Median = EazyMZDataProcess.FD_Line(ABS_FD_List)
                    FD_Median_N = FD_Median*(-1)
                    Diff_Median = EazyMZDataProcess.Diff_Line(Diff_List)
                    SD_Median = EazyMZDataProcess.SD_Line(SD_List)
                    FD_P = list(filter(lambda x: FD_List[x] > FD_Median, range(1, len(FD_List))))
                    FD_N = list(filter(lambda x: FD_List[x] < FD_Median_N, range(1, len(FD_List))))
                    FD_Change = list(filter(lambda x: FD_List[x-1] > 0 and FD_List[x] < 0, range(1, len(FD_List))))
                    Diff_P = list(filter(lambda x: Diff_List[x] > Diff_Median, range(len(Diff_List))))
                    Diff_N = list(filter(lambda x: Diff_List[x] < (Diff_Median*-1), range(len(Diff_List))))
                    SD_Place = list(filter(lambda x: SD_List[x] < SD_Median and SD_List[x] < 0 and (SD_List[x-1]>SD_List[x] or SD_List[x+1]>SD_List[x]), range(1,len(SD_List)-1)))
                    if len(Diff_P) == 0 or len(FD_P) == 0 or len(Diff_N) == 0 or len(FD_N) == 0 or len(SD_Place) == 0:
                        return []
                    else:
                        # 引入人机交互，调整参数？
                        Begin_Place,Complet_BP = EazyMZDataProcess.Find_FContinuous(FD_List, FD_P, Diff_P, FDMode='P',mergeRule=self.get_param('MergeRule'),FC_Number=self.get_param('FeatureDetectPlot'))
                        Begin_Place += 2
                        Complet_BP += 2
                        End_Place,Complet_EP = EazyMZDataProcess.Find_FContinuous(FD_List, FD_N, Diff_N, FDMode='N',mergeRule=self.get_param('MergeRule'),FC_Number=self.get_param('FeatureDetectPlot'))
                        End_Place += 2
                        Complet_EP += 2
                        Peak_Place = EazyMZDataProcess.Find_SDChange(FD_Change, SD_Place)  # 峰顶位置
                        Peak_Place += 2
                        #Peak_Place = np.array(list(filter(lambda x:Final_RTL <= Auto_RT_List[x] <= Final_RTR,Peak_Place)))
                    if len(Begin_Place) >= 1 and len(End_Place) >= 1 and Begin_Place[0] < End_Place[-1]:
                        FD_Sequence = np.concatenate((Begin_Place, End_Place))
                        FD_Sequence.sort()
                        Peak_Begin = [Begin_Place[0]]
                        Peak_End = []
                        for i_FD_Seq in range(1, len(Begin_Place)):
                            Begin_Value = Begin_Place[i_FD_Seq]
                            temp_BValue_Place = np.where(FD_Sequence == Begin_Value)[0]
                            temp_EValue_Place = np.where(FD_Sequence == Peak_Begin[-1])[0]
                            if FD_Sequence[temp_BValue_Place[0]] > Peak_Begin[-1] and FD_Sequence[temp_BValue_Place[0]-1] in Begin_Place.tolist():
                                Peak_Begin[-1] = FD_Sequence[temp_BValue_Place[0]]
                            elif temp_EValue_Place[0] < len(FD_Sequence)-1 and temp_BValue_Place[0] < len(FD_Sequence)-1:
                                if FD_Sequence[temp_EValue_Place[0]+1] > Peak_Begin[-1] and FD_Sequence[temp_EValue_Place[0]+1] in End_Place.tolist():
                                    Peak_End.append(FD_Sequence[temp_EValue_Place[0]+1])
                                    Peak_Begin.append(FD_Sequence[temp_BValue_Place[0]])                           
                        if len(Peak_End) > 0:
                            temp_End_List = list(filter(lambda x: x > Peak_End[-1] and x > Peak_Begin[-1], End_Place))
                            if len(temp_End_List) > 0:
                                temp_End_Place = EazyMZDataProcess.ClosestPosition(Peak_Begin[-1], temp_End_List)
                                Peak_End.append(temp_End_List[temp_End_Place])
                            else:
                                del Peak_Begin[-1]
                        elif len(Peak_End) == 0:
                            temp_End_List = list(filter(lambda x: x > Peak_Begin[-1], End_Place))
                            if len(temp_End_List) > 0:
                                temp_End_Place = EazyMZDataProcess.ClosestPosition(Peak_Begin[-1], temp_End_List)
                                Peak_End.append(temp_End_List[temp_End_Place])
                            else:
                                del Peak_Begin[-1]
                        # 删除上升和下降之间差距过大的峰，例如冲顶峰和平峰以及噪声
                        gap_del = []
                        for i_com in range(len(Peak_Begin)):
                            BP_gap = np.where(Begin_Place==Peak_Begin[i_com])[0][0]
                            EP_gap = np.where(End_Place==Peak_End[i_com])[0][0]
                            if Complet_EP[EP_gap]-Complet_BP[BP_gap]>=self.get_param('UpDown_gap'):
                                gap_del.append(i_com)
                        Peak_Begin = np.array(Peak_Begin)
                        Peak_Begin = list(np.delete(Peak_Begin,gap_del))
                        Peak_End = np.array(Peak_End)
                        Peak_End = list(np.delete(Peak_End,gap_del))
                        if len(Peak_Begin) != len(Peak_End):
                            raise ValueError('Peak len not match')
                        Peak_Top = []
                        for i_top in range(len(Peak_Begin)-1, -1, -1):
                            temp_top = np.where((Peak_Begin[i_top] < Peak_Place) & (Peak_End[i_top] > Peak_Place))[0]
                            if len(temp_top) >= 1:
                                top_ran = np.array(range(len(temp_top)))
                                Auto_Int_List = np.array(Auto_Int_List)
                                top_ran = list(filter(lambda x: Auto_Int_List[Peak_Place[temp_top[x]]] == np.max(Auto_Int_List[Peak_Place[temp_top]]), top_ran))
                                Peak_Top.append(Peak_Place[temp_top[top_ran[0]]])
                            else:
                                del Peak_Begin[i_top]
                                del Peak_End[i_top]
                        Peak_Top = list(map(lambda x: Peak_Top[-x], range(1, len(Peak_Top)+1)))
                        unfit = []
                        for i_top in range(len(Peak_Begin)):
                            if Peak_Begin[i_top]!=0:
                                if Peak_Begin[i_top] < self.get_param('Flow_RT'):
                                    begin_ran = np.array(range(self.get_param('Flow_RT')))
                                else:
                                    begin_ran = np.array(range(Peak_Begin[i_top]-self.get_param('Flow_RT')+1, Peak_Begin[i_top]+1))
                                temp_be = list(filter(lambda x: Auto_Int_List[x] == np.min(Auto_Int_List[begin_ran]), begin_ran))
                                Peak_Begin[i_top] = temp_be[-1]
                            if Peak_End[i_top]!=len(Auto_Int_List)-1:
                                if Peak_End[i_top]+self.get_param('Flow_RT') >= len(Auto_Int_List):
                                    end_ran = np.array(range(len(Auto_Int_List)-self.get_param('Flow_RT'), len(Auto_Int_List)))
                                else:
                                    end_ran = np.array(range(Peak_End[i_top], Peak_End[i_top]+self.get_param('Flow_RT')))
                                temp_be = list(filter(lambda x: Auto_Int_List[x] == np.min(Auto_Int_List[end_ran]), end_ran))
                                Peak_End[i_top] = temp_be[0]
                            top_ran = np.array(range(Peak_Top[i_top]-1, Peak_Top[i_top]+2))
                            temp_be = list(filter(lambda x: Auto_Int_List[x] == np.max(Auto_Int_List[top_ran]), top_ran))
                            Peak_Top[i_top] = temp_be[0]
                            if Peak_Begin[i_top]>=Peak_End[i_top]:
                                unfit.append(i_top)
                        if len(unfit)>0:
                            Peak_Begin = np.array(Peak_Begin)
                            Peak_End = np.array(Peak_End)
                            Peak_Top = np.array(Peak_Top)
                            Peak_Begin = list(np.delete(Peak_Begin,unfit))
                            Peak_End = list(np.delete(Peak_End,unfit))
                            Peak_Top = list(np.delete(Peak_Top,unfit))
                        return Peak_Begin,Peak_End,Peak_Top,Origin_Int_List,Auto_RT_List,MZ
                    else:
                        return []
                else:
                    return []
            else:
                return []
        else:
            return []
        
    def FPD_importer(path):
        def str_num(x):
            x = re.split(' |\\|[|]',x)
            x = list(filter(lambda x:len(x)>0,x))
            x = list(map(lambda x:float(x),x))
            return x
        rawData = pd.read_excel(path)
        rawData['RTList'] = rawData['RTList'].apply(lambda x:str_num(x[1:len(x)-1]))
        rawData['IntList'] = rawData['IntList'].apply(lambda x:str_num(x[1:len(x)-1]))
        return rawData
    def RI_to_RT(RI,RIIS):
        n_place = list(filter(lambda x:RIIS.at[x,'C']<=RI/100,range(len(RIIS))))
        n1_place = list(filter(lambda x:RIIS.at[x,'C']>RI/100,range(len(RIIS))))
        if len(n_place)==0:
            n_place = n1_place[0]
            n1_place = n1_place[1]
        elif len(n1_place)==0:
            n1_place = n_place[-1]
            n_place = n_place[-2] 
        else:
            n_place = n_place[-1]
            n1_place = n1_place[0]
        C_number = RIIS['C'][n_place]
        t_n = RIIS['RT'][n_place]
        t_n1 = RIIS['RT'][n1_place]
        RT = (RI/100-C_number)/(RIIS['C'][n1_place]-RIIS['C'][n_place])*(t_n1-t_n)+t_n
        return RT
    def RT_to_RI(RT,RIIS):
        n_place = list(filter(lambda x:RIIS.at[x,'RT']<=RT,range(len(RIIS))))
        n1_place = list(filter(lambda x:RIIS.at[x,'RT']>RT,range(len(RIIS))))
        if len(n_place)==0:
            n_place = n1_place[0]
            n1_place = n1_place[1]
        elif len(n1_place)==0:
            n1_place = n_place[-1]
            n_place = n_place[-2] 
        else:
            n_place = n_place[-1]
            n1_place = n1_place[0]
        C_number = RIIS['C'][n_place]
        t_n = RIIS['RT'][n_place]
        t_n1 = RIIS['RT'][n1_place]
        RI = 100*(C_number+(RIIS['C'][n1_place]-RIIS['C'][n_place])*(RT-t_n)/(t_n1-t_n))
        return RI
    def Heuristic_PeakDetect(self):
        warnings.filterwarnings('ignore')
        def RI_to_RT(RI,RIIS):
            n_place = list(filter(lambda x:RIIS.at[x,'C']<=RI/100,range(len(RIIS))))
            n1_place = list(filter(lambda x:RIIS.at[x,'C']>RI/100,range(len(RIIS))))
            if len(n_place)==0:
                n_place = n1_place[0]
                n1_place = n1_place[1]
            elif len(n1_place)==0:
                n1_place = n_place[-1]
                n_place = n_place[-2] 
            else:
                n_place = n_place[-1]
                n1_place = n1_place[0]
            C_number = RIIS['C'][n_place]
            t_n = RIIS['RT'][n_place]
            t_n1 = RIIS['RT'][n1_place]
            RT = (RI/100-C_number)/(RIIS['C'][n1_place]-RIIS['C'][n_place])*(t_n1-t_n)+t_n
            return RT
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        def gaussian(x,a=1,b=10,c=2):
            y=a*math.exp((-1*(x-b)**2)/(2*c**2))
            return y
        def RT_to_RI(RT,RIIS):
            n_place = list(filter(lambda x:RIIS.at[x,'RT']<=RT,range(len(RIIS))))
            n1_place = list(filter(lambda x:RIIS.at[x,'RT']>RT,range(len(RIIS))))
            if len(n_place)==0:
                n_place = n1_place[0]
                n1_place = n1_place[1]
            elif len(n1_place)==0:
                n1_place = n_place[-1]
                n_place = n_place[-2] 
            else:
                n_place = n_place[-1]
                n1_place = n1_place[0]
            C_number = RIIS['C'][n_place]
            t_n = RIIS['RT'][n_place]
            t_n1 = RIIS['RT'][n1_place]
            RI = 100*(C_number+(RIIS['C'][n1_place]-RIIS['C'][n_place])*(RT-t_n)/(t_n1-t_n))
            return RI
        gaussian_params = [0.5, 1, 2, 3, 4, 5]
        Heuristic_peak_list = self.Heuristic_peak_list.copy()
        Final_Peak_Detect = pd.DataFrame(columns=['AverageMZ', 'RT', 'RI', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
        Heuristic_peak_list['ActualNumber'] = 0
        for ii in range(len(Heuristic_peak_list)):   
            MZ = Heuristic_peak_list.at[ii,'AverageMZ']
            RI_left = Heuristic_peak_list.at[ii,'RI_left']
            RI_right = Heuristic_peak_list.at[ii,'RI_right']
            RI_List = Heuristic_peak_list.at[ii,'RI_List']
            RTL = RI_to_RT(RI_left,self.RIIS)
            RTR = RI_to_RT(RI_right,self.RIIS)
            #Gradient_Code = Heuristic_peak_list.at[ii,'Gradient Code']
            [RT_List,Int_List_Origin] = self.ExtractDataPoint(MZ,RTR+self.get_param('RT_Tor'),RTL-self.get_param('RT_Tor'),smooth_index=0)
            if len(RT_List) > self.get_param('Points'):
                Peak_Result = self.PeakDetecter_for_Heuristic(RT_List,Int_List_Origin,MZ,WhetherAppend=True)
                if len(Peak_Result) >0:
                    Peak_Begin = Peak_Result[0]
                    Peak_End = Peak_Result[1]
                    Peak_Top = Peak_Result[2]
                    Auto_Int_List = Peak_Result[3]
                    Auto_RT_List = Peak_Result[4]
                    MZ = Peak_Result[5]
                    results = []
                    RT_match = []
                    for v in range(len(Peak_Begin)):
                        #temp_Final_Peak = pd.DataFrame(columns=['AverageMZ', 'RT', 'RI', 'Int','RTList','SimilarityScore'])
                        RT = Auto_RT_List[Peak_Top[v]]
                        RI = RT_to_RI(RT, self.RIIS)
                        IntList = Auto_Int_List[Peak_Begin[v]:Peak_End[v]+1]
                        RTList = Auto_RT_List[Peak_Begin[v]:Peak_End[v]+1]
                        Int = max(IntList)
                        if IntList[0]/Int <= 0.20 or IntList[-1]/Int <= 0.20:
                            if RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points') and len(list(filter(lambda x:(1-self.get_param('RI_Tor'))<=RI/x<=(1+self.get_param('RI_Tor')),RI_List)))>=1: 
                                Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                SimilarityScore = max(SimilarityScores)
                                temp_Final_Peak = pd.DataFrame([(MZ, RT, RI, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'RI', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                results.append(temp_Final_Peak)
                                RT_match.append(RT)
                            elif RT < self.RIIS.at[0,'RT'] and RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points') and len(list(filter(lambda x:(1-self.get_param('RI_Tor')*2)<=RI/x<=(1+self.get_param('RI_Tor')*2),RI_List)))>=1:
                                Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                SimilarityScore = max(SimilarityScores)
                                temp_Final_Peak = pd.DataFrame([(MZ, RT, RI, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'RI', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                results.append(temp_Final_Peak)
                                RT_match.append(RT)
                    if len(results) < Heuristic_peak_list.at[ii,'PeakNumber'] and len(Peak_Begin) > 0:
                        Peak_End.insert(0,0)
                        Peak_Begin.append(len(Auto_RT_List))
                        for v in range(len(Peak_Begin)):
                            if Peak_Begin[v]-Peak_End[v] >= self.get_param('Points') and max(Auto_Int_List[Peak_End[v]:Peak_Begin[v]]) >= self.get_param('min_Int'):
                                Peak_Result = self.PeakDetecter_for_Heuristic(Auto_RT_List[Peak_End[v]:Peak_Begin[v]],Auto_Int_List[Peak_End[v]:Peak_Begin[v]],MZ,WhetherAppend=False)
                                if len(Peak_Result) >0:
                                    Peak_Begin_Iterate = Peak_Result[0]
                                    Peak_End_Iterate = Peak_Result[1]
                                    Peak_Top_Iterate = Peak_Result[2]
                                    Auto_Int_List_Iterate = Peak_Result[3]
                                    Auto_RT_List_Iterate = Peak_Result[4]
                                    MZ_Iterate = Peak_Result[5]
                                    for vi in range(len(Peak_Begin_Iterate)):
                                        RT = Auto_RT_List_Iterate[Peak_Top_Iterate[vi]]
                                        RI = RT_to_RI(RT, self.RIIS)
                                        IntList = Auto_Int_List_Iterate[Peak_Begin_Iterate[vi]:Peak_End_Iterate[vi]+1]
                                        RTList = Auto_RT_List_Iterate[Peak_Begin_Iterate[vi]:Peak_End_Iterate[vi]+1]
                                        Int = max(IntList)
                                        if IntList[0]/Int <= 0.20 or IntList[-1]/Int <= 0.20 and len(list(filter(lambda x:abs(RT-x)<self.get_param('RT_Tor')/3,RT_match)))<0:
                                            if RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points') and len(list(filter(lambda x:(1-self.get_param('RI_Tor'))<=RI/x<=(1+self.get_param('RI_Tor')),RI_List)))>=1: 
                                                Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                                SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                                SimilarityScore = max(SimilarityScores)
                                                temp_Final_Peak = pd.DataFrame([(MZ_Iterate, RT, RI, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'RI', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                                results.append(temp_Final_Peak)
                                                RT_match.append(RT)
                                            elif RT < self.RIIS.at[0,'RT'] and RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points') and len(list(filter(lambda x:(1-self.get_param('RI_Tor')*2)<=RI/x<=(1+self.get_param('RI_Tor')*2),RI_List)))>=1:
                                                Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                                SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                                SimilarityScore = max(SimilarityScores)
                                                temp_Final_Peak = pd.DataFrame([(MZ_Iterate, RT, RI, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'RI', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                                results.append(temp_Final_Peak)
                                                RT_match.append(RT)
                    if len(results)>0:
                        temp_Final_Peak = pd.concat(results)
                        temp_Final_Peak.reset_index(drop=True,inplace=True)
                        Int_max = max(temp_Final_Peak['Int'])
                        Score_matrix= np.zeros([len(RI_List),len(temp_Final_Peak)])
                        for i_Sm_row in range(len(RI_List)):
                            for i_Sm_col in range(len(temp_Final_Peak)):
                                if abs(temp_Final_Peak.at[i_Sm_col,'RI']/RI_List[i_Sm_row]-1) <= self.get_param('RI_Tor'):
                                    Score_matrix[i_Sm_row,i_Sm_col] = 1-abs(temp_Final_Peak.at[i_Sm_col,'RI']/RI_List[i_Sm_row]-1)/self.get_param('RI_Tor')+temp_Final_Peak.at[i_Sm_col,'Int']/Int_max
                        Sm_row,Sm_col = linear_sum_assignment(Score_matrix,True)
                        Valified_TFP = []
                        for i_Sm_row,i_Sm_col in zip(Sm_row,Sm_col):
                            if Score_matrix[i_Sm_row,i_Sm_col] != 0:
                                Valified_TFP.append(i_Sm_col)
                        temp_Final_Peak = temp_Final_Peak.iloc[Valified_TFP,:].copy()
                        temp_Final_Peak.reset_index(drop=True,inplace=True)
                        if len(temp_Final_Peak)<= Heuristic_peak_list.at[ii,'PeakNumber']:
                            Heuristic_peak_list.loc[ii,'ActualNumber'] = len(temp_Final_Peak)
                        else:
                            Heuristic_peak_list.loc[ii,'ActualNumber'] = Heuristic_peak_list.at[ii,'PeakNumber'].copy()
                            temp_Final_Peak.sort_values('Int',inplace=True,ignore_index=True)
                            temp_Final_Peak = temp_Final_Peak.iloc[range(Heuristic_peak_list.at[ii,'PeakNumber']),:].copy()
                        Final_Peak_Detect = pd.concat([Final_Peak_Detect, temp_Final_Peak],ignore_index=True)
        data_points = np.column_stack((Final_Peak_Detect['RT'],Final_Peak_Detect['AverageMZ']))
        radius = self.get_param('RT_Tor')/3
        kdtree = KDTree(data_points)
        del_list = []
        for ii in range(len(Final_Peak_Detect)):   
            Match_List = kdtree.query_ball_point([Final_Peak_Detect.at[ii,'RT'],Final_Peak_Detect.at[ii,'AverageMZ']], r=radius)
            Match_List = list(filter(lambda x:x > ii and abs(Final_Peak_Detect.at[ii,'AverageMZ']-Final_Peak_Detect.at[x,'AverageMZ'])/Final_Peak_Detect.at[ii,'AverageMZ']< self.get_param('MS1_Tor'),Match_List))
            for iii in Match_List:
                if Final_Peak_Detect.at[ii,'Int'] >= Final_Peak_Detect.at[iii,'Int']:
                    del_list.append(iii)
                else:
                    del_list.append(ii)
        del_list = list(set(del_list))
        Final_Peak_Detect.drop(del_list,inplace=True)
        Final_Peak_Detect.reset_index(drop=True,inplace=True)
        self.Final_Peak_Detect = Final_Peak_Detect.copy()     
        return Final_Peak_Detect.copy()
    
    def set_param(self,name,value):
        if name in self.__param:
            self.__param[name]=value
        else:
            print('no such param')
    def get_param(self,key=''):
        if len(key) == 0:
            print(self.__param)
        else:
            return self.__param[key]
    def get_MZ(self,x):
        return self.Final_Peak_Detect.at[x,'AverageMZ']
    def get_RT(self,x):
        return self.Final_Peak_Detect.at[x,'RT']
    def TFFD(x, Auto_Int_List, Auto_RT_List):
        FirstDerivative = (Auto_Int_List[x+1]*8+Auto_Int_List[x-2]-Auto_Int_List[x-1] * 8-Auto_Int_List[x+2])/((Auto_RT_List[x+2]-Auto_RT_List[x-2])*3)
        return FirstDerivative
    def TFSD(x, Auto_Int_List, Auto_RT_List):
        SecondDerivative = (Auto_Int_List[x+1]*16+Auto_Int_List[x-1]*16-Auto_Int_List[x] * 30-Auto_Int_List[x+2]-Auto_Int_List[x-2])/((Auto_RT_List[x+2]-Auto_RT_List[x-2])*3)
        return SecondDerivative
    def FD_Line(ABS_FD_List):
        ABS_FD_List = list(filter(lambda x: x < max(ABS_FD_List)*0.20 and x > 0, ABS_FD_List))
        FD_Median = np.median(ABS_FD_List)
        return FD_Median
    def SD_Line(SD_List):
        #SD_List = abs(SD_List)
        if len(SD_List) >= 1:
            SD_Median = np.median(SD_List[SD_List<0])
            return SD_Median
        else:
            return 0
    def Diff_Line(Diff_List,limit=0.1):
        Diff_List = list(map(lambda x: abs(x), Diff_List))
        Diff_List = list(filter(lambda x: x < max(Diff_List)*limit and x != 0, Diff_List))
        if len(Diff_List) >= 1:
            Diff_Median = np.median(Diff_List)
        else:
            Diff_Median = 0
        return Diff_Median 
    def Find_FContinuous(FD_List,FD_P, Diff_P, FDMode='P' ,FC_Number=2 ,mergeRule='Intersection'):
        if mergeRule == 'Intersection':
            First_List = list(filter(lambda x:x in FD_P,Diff_P))
        elif mergeRule == 'Union':
            First_List = list(set(FD_P+Diff_P))
            First_List.sort()
        if len(First_List)>0:
            if FDMode == 'P':
                Second_List = [First_List[0]]
                Second_Score = [First_List[0]]
                for i_FSB in range(1,len(First_List)):
                    if First_List[i_FSB]-First_List[i_FSB-1]<=FC_Number:
                        Second_Score[-1] = First_List[i_FSB]
                    else:
                        Second_List.append(First_List[i_FSB])
                        Second_Score.append(First_List[i_FSB])
                Second_List = np.array(Second_List)
                Second_Score = np.array(Second_Score)
                temp_Place = np.where(Second_Score-Second_List >= FC_Number*2-1)[0]
            elif FDMode == 'N':
                Second_List = [First_List[0]]
                Second_Score = [First_List[0]]
                for i_FSB in range(1,len(First_List)):
                    if First_List[i_FSB]-First_List[i_FSB-1]<=FC_Number:
                        Second_Score[-1] = First_List[i_FSB]
                    else:
                        Second_List.append(First_List[i_FSB])
                        Second_Score.append(First_List[i_FSB])
                Second_List = np.array(Second_List)
                Second_Score = np.array(Second_Score)
                temp_Place = np.where(Second_Score-Second_List >= FC_Number*2-1)[0]
                Second_List = Second_List[temp_Place]
                Second_Score = Second_Score[temp_Place]
                for i_FSB in range(len(Second_List)):
                    FD_List_short = FD_List[Second_List[i_FSB]:Second_Score[i_FSB]]
                    Int_min = min(FD_List_short)
                    count = 0
                    for ii_FSB in range(np.where(FD_List_short==Int_min)[0][0],len(FD_List_short)):
                        if FD_List_short[ii_FSB] >= Int_min * 0.05:
                            count += 1
                        if count >= 3 and ii_FSB >= FC_Number*2-1:
                            Second_Score[i_FSB] = Second_List[i_FSB]+ii_FSB
                            break
            if len(temp_Place) > 0:
                if FDMode == 'P':
                    return Second_List[temp_Place],Second_Score[temp_Place]
                elif FDMode == 'N':
                    return Second_Score,Second_List
            else:
                return np.array([]),np.array([])
        else:
            return np.array([]),np.array([])
    def Find_SDChange(FD_Change, SD_Place):
        if len(FD_Change) <= len(SD_Place):
            list_a = np.array(FD_Change)
            list_b = np.array(SD_Place)
        else:
            list_b = np.array(FD_Change)
            list_a = np.array(SD_Place)
        Change_Place = []
        for i in range(len(list_a)):
            temp = np.where(abs(list_a[i]-list_b) <= 2)[0]
            if len(temp) >= 1:
                Change_Place.append(i)
        if len(Change_Place) > 0:
            Change_Place = np.array(Change_Place)
            return list_a[Change_Place]
        else:
            return np.array([])
    def Baseline(self,MZ, RTR, RTL, s_min=False):
        if s_min == True:
            RTR = RTR * 60          # min -> s
            RTL = RTL * 60
        [RT_List_L,Int_List_L] = EazyMZDataProcess.ExtractDataPoint(self,MZ,RTL,2*RTL-1*RTR,smooth_index=0)
        [RT_List_R,Int_List_R] = EazyMZDataProcess.ExtractDataPoint(self,MZ,2*RTR-1*RTL,RTR,smooth_index=0)
        if len(Int_List_R)>0 and len(Int_List_L)>0:
            BaselineValue = min(0.5*(np.median(Int_List_R)+sum(Int_List_R)/len(Int_List_R)),0.5*(np.median(Int_List_L)+sum(Int_List_L)/len(Int_List_L)))
        elif len(Int_List_R)>0 and len(Int_List_L)==0:
            BaselineValue = 0.5*(np.median(Int_List_R)+sum(Int_List_R)/len(Int_List_R))
        elif len(Int_List_R)==0 and len(Int_List_L)>0:
            BaselineValue = 0.5*(np.median(Int_List_L)+sum(Int_List_L)/len(Int_List_L))
        else:
            BaselineValue = 1
        Int_List_L_diff = list(map(lambda x:x-BaselineValue,Int_List_L))
        Int_List_L_diff = list(filter(lambda x:x>0,Int_List_L_diff))
        Int_List_R_diff = list(map(lambda x:x-BaselineValue,Int_List_R))
        Int_List_R_diff = list(filter(lambda x:x>0,Int_List_R_diff))
        if len(Int_List_L_diff)>0 and len(Int_List_R_diff)>0:
            Noise = min(max(np.median(Int_List_R_diff),sum(Int_List_R_diff)/len(Int_List_R_diff)),max(np.median(Int_List_L_diff),sum(Int_List_L_diff)/len(Int_List_L_diff)))
            if Noise<=0:
                Noise = 1
        else:
            Noise = 1
        return BaselineValue,Noise
    def Calculate_SN(self,drop=False,Threshold=10):
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        bar = Bar('Calculate_SN ',max=len(self.Final_Peak_Detect))
        SN = []
        for i in range(len(self.Final_Peak_Detect)):
            bar.next()
            [RT_List,Int_List] = self.ExtractDataPoint(self.Final_Peak_Detect.at[i,'AverageMZ'],self.Origin_RT_List[-1],0,smooth_index=0)
            Int_List = np.array(Int_List)
            Auto_Int_List = np.array(list(map(lambda x:GaussSmooth(Int_List[x-2:x+3]) if 2<=x<=len(Int_List)-2 else Int_List[x],range(len(Int_List)))))
            if np.log10(max(Int_List)/self.get_param('min_Int'))>=2:
                noise = sum(abs(Int_List[(Int_List>0)&(Int_List<max(Int_List)*0.1)]-Auto_Int_List[(Int_List>0)&(Int_List<max(Int_List)*0.1)]))/len(Int_List[(Int_List>0)&(Int_List<max(Int_List)*0.1)])
            else:
                noise = sum(abs(Int_List[Int_List>0]-Auto_Int_List[Int_List>0]))/len(Int_List[Int_List>0])
            SN.append(self.Final_Peak_Detect.at[i,'Int']/noise)
        bar.finish()
        self.Final_Peak_Detect['S/N'] = SN
        if drop == True:
            self.Final_Peak_Detect = self.Final_Peak_Detect[self.Final_Peak_Detect['S/N']>Threshold].copy()
            self.Final_Peak_Detect.reset_index(drop=True,inplace=True)
    def set_RIIS(self,RIIS):
        self.RIIS = pd.read_excel(RIIS)
        if 'RT' not in list(self.RIIS.keys()):
            self.RIIS['RT']=0.00
            self.RIIS['Int']=0.00
            self.RIIS['Candidate_RT'] = ''
            self.RIIS['Candidate_Int'] = ''
            RTL = self.get_param('min_RT')
            RTR = self.Origin_RT_List[-1]
            for i in range(len(self.RIIS)):                
                MZ = self.RIIS.at[i,'m/z']
                [RT_List, Int_List]=self.ExtractDataPoint(MZ,RTR,RTL)
                RT_place_1 = Int_List.index(max(Int_List))
                Int_List_d = Int_List[0:RT_place_1-int(self.get_param('Points'))*2]+Int_List[RT_place_1+int(self.get_param('Points'))*2:-1]
                RT_place_2 = Int_List_d.index(max(Int_List_d))
                if RT_place_2 >= RT_place_1-int(self.get_param('Points')):
                    RT_place_2 += int(self.get_param('Points'))*4
                if Int_List[RT_place_1]>=self.__param['min_Int']:
                    self.RIIS.at[i,'Candidate_RT'] = [RT_List[RT_place_1],RT_List[RT_place_2]]
                    self.RIIS.at[i,'Candidate_Int'] = [Int_List[RT_place_1],Int_List[RT_place_2]]
                    self.RIIS.at[i,'RT'] = RT_List[RT_place_1].copy()
                    self.RIIS.at[i,'Int'] = Int_List[RT_place_1].copy()
                if Int_List[RT_place_1]/Int_List[RT_place_2] >=5:
                    RTL = RT_List[RT_place_1]
            del_list = list(filter(lambda x:self.RIIS.at[x,'RT']==0,range(len(self.RIIS))))
            self.RIIS.drop(del_list,inplace=True)
            self.RIIS.reset_index(drop=True,inplace=True)   
            del_fix = 0
            for i in range(len(self.RIIS.iloc[:,0])): #len(self.RIIS.iloc[:,0])
                i = i - del_fix
                if i == 0:
                    if self.RIIS.at[i,'RT'] >= self.RIIS['RT'][i+1] and len(self.RIIS['Candidate_RT'][i]) == 2:
                        self.RIIS.at[i,'RT'] = self.RIIS.at[i,'Candidate_RT'][1]
                        if self.RIIS.at[i,'RT'] >= self.RIIS.at[i+1,'RT']:
                            self.RIIS.drop(i,inplace=True)
                            self.RIIS.reset_index(drop=True,inplace=True)  
                            del_fix = del_fix + 1
                    elif self.RIIS.at[i,'RT'] >= self.RIIS.at[i+1,'RT'] and len(self.RIIS.at[i,'Candidate_RT']) == 1:
                        self.RIIS.drop(i,inplace=True)
                        self.RIIS.reset_index(drop=True,inplace=True)  
                        del_fix = del_fix + 1
                elif i == len(self.RIIS.iloc[:,0])-1:
                    if self.RIIS.at[i,'RT'] <= self.RIIS.at[i-1,'RT'] and len(self.RIIS.at[i,'Candidate_RT']) == 2:
                        self.RIIS.at[i,'RT'] = self.RIIS.at[i,'Candidate_RT'][1]
                        if self.RIIS.at[i,'RT'] <= self.RIIS.at[i-1,'RT']:
                            self.RIIS.drop(i,inplace=True)
                            self.RIIS.reset_index(drop=True,inplace=True)  
                            del_fix = del_fix + 1
                    elif self.RIIS.at[i,'RT'] <= self.RIIS.at[i-1,'RT'] and len(self.RIIS.at[i,'Candidate_RT']) == 1:
                        self.RIIS.drop(i,inplace=True)
                        self.RIIS.reset_index(drop=True,inplace=True)  
                        del_fix = del_fix + 1
                else:
                    if self.RIIS.at[i,'RT'] > self.RIIS.at[i-1,'RT'] and self.RIIS.at[i,'RT'] < self.RIIS.at[i+1,'RT']:
                        continue
                    if self.RIIS.at[i,'RT'] < self.RIIS.at[i-1,'RT']:
                        if len(self.RIIS.at[i,'Candidate_RT']) == 2:
                            self.RIIS.at[i,'RT'] = self.RIIS.at[i,'Candidate_RT'][1]
                            if self.RIIS.at[i,'RT'] <= self.RIIS.at[i-1,'RT'] or self.RIIS.at[i,'RT'] >=self.RIIS.at[i+1,'RT']:
                                self.RIIS.drop(i,inplace=True)
                                self.RIIS.reset_index(drop=True,inplace=True)  
                                del_fix = del_fix + 1 
                    elif self.RIIS.at[i,'RT'] > self.RIIS.at[i+1,'RT']:
                        if i >=2 :
                            i_n = ((self.RIIS.at[i,'RT']-self.RIIS.at[i-1,'RT'])/(self.RIIS.at[i,'C']-self.RIIS.at[i-1,'C']))/((self.RIIS.at[i-1,'RT']-self.RIIS.at[i-2,'RT'])/(self.RIIS.at[i-1,'C']-self.RIIS.at[i-2,'C']))
                            i_n_1 = ((self.RIIS.at[i+1,'RT']-self.RIIS.at[i-1,'RT'])/(self.RIIS.at[i+1,'C']-self.RIIS.at[i-1,'C']))/((self.RIIS.at[i-1,'RT']-self.RIIS['RT'][i-2])/(self.RIIS['C'][i-1]-self.RIIS['C'][i-2]))
                            if abs(1-i_n) > abs(1-i_n_1):
                                if len(self.RIIS['Candidate_RT'][i]) == 2:
                                        self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][1]
                                        if self.RIIS['RT'][i] <= self.RIIS['RT'][i-1] or self.RIIS['RT'][i] >=self.RIIS['RT'][i+1]:
                                            self.RIIS.drop(i,inplace=True)
                                            self.RIIS.reset_index(drop=True,inplace=True)  
                                            del_fix = del_fix + 1 
                        elif len(self.RIIS['Candidate_RT'][i]) == 2:
                                self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][1]
                                if self.RIIS['RT'][i] <= self.RIIS['RT'][i-1] or self.RIIS['RT'][i] >=self.RIIS['RT'][i+1]:
                                    self.RIIS.drop(i,inplace=True)
                                    self.RIIS.reset_index(drop=True,inplace=True)  
                                    del_fix = del_fix + 1      
            for i in range(1,len(self.RIIS.iloc[:,1])-1):
                if len(self.RIIS['Candidate_RT'][i]) == 2:
                    if min(self.RIIS['Candidate_Int'][i])/max(self.RIIS['Candidate_Int'][i])>0.4 and self.RIIS['RT'][i] != self.RIIS['Candidate_RT'][i][1]:
                        if self.RIIS.at[i-1,'RT']<self.RIIS['Candidate_RT'][i][1]<self.RIIS.at[i+1,'RT']:
                            score_1 = (self.RIIS['Candidate_RT'][i][0]-self.RIIS.at[i-1,'RT'])/(self.RIIS.at[i+1,'RT']-self.RIIS['Candidate_RT'][i][0])
                            score_2 = (self.RIIS['Candidate_RT'][i][1]-self.RIIS.at[i-1,'RT'])/(self.RIIS.at[i+1,'RT']-self.RIIS['Candidate_RT'][i][1])
                            if abs(1-score_2)<abs(1-score_1):
                                self.RIIS.at[i,'RT']=self.RIIS.at[i,'Candidate_RT'][1]
        else:
            if max(self.RIIS['RT'])<100:
                self.RIIS['RT'] = self.RIIS['RT'].apply(lambda x:x*60)
            for i in range(len(self.RIIS)): 
                MZ = self.RIIS.at[i,'m/z']
                [RT_List, Int_List]=self.ExtractDataPoint(MZ,self.RIIS.at[i,'RT']+self.get_param('RT_Tor'),self.RIIS.at[i,'RT']-self.get_param('RT_Tor'),smooth_index=self.get_param('smooth'))
                RT_place = Int_List.index(max(Int_List))
                if max(Int_List) > self.get_param('min_Int'):
                    RT_place = Int_List.index(max(Int_List))
                    self.RIIS.at[i,'RT'] = RT_List[RT_place]
                else:
                    continue
                
    def Calculate_RI(self):
        if 'RT' not in self.RIIS.keys():
            self.RIIS['Candidate_RT'] = ''
            self.RIIS['Candidate_Int'] = ''
            self.RIIS['RT']=np.nan
            Missing_RI = []
            for i in range(len(self.RIIS.iloc[:,1])):
                IS_place = np.where(abs(self.Final_Peak_Detect['AverageMZ']-self.RIIS.iloc[i,1])/self.RIIS.iloc[i,1]<self.__param['MS1_Tor'])[0]
                Candidate_RT = []
                Candidate_Int = []
                if len(IS_place)>0:
                    for ii in IS_place:
                        Candidate_RT.append(self.Final_Peak_Detect['RT'][ii])
                        Candidate_Int.append(self.Final_Peak_Detect['Int'][ii])
                    temp_DF = pd.DataFrame({'RT':Candidate_RT,'Int':Candidate_Int})
                    temp_DF.sort_values(by='Int',ascending=False,axis=0,inplace=True)
                    temp_DF.reset_index(drop=True,inplace=True)
                    self.RIIS['Candidate_RT'][i] = list(temp_DF['RT'])
                    self.RIIS['Candidate_Int'][i] = list(temp_DF['Int'])
                    self.RIIS['RT'][i] = temp_DF['RT'][0]
                else:
                    self.RIIS['RT'][i] = np.nan
                    Missing_RI.append(i)
            if len(np.where(np.isnan(self.RIIS['RT'])==True)[0])/len(self.RIIS['RT'])>0.4:
                raise ValueError('Too much missing RI stander')
            else:
                self.RIIS.drop(Missing_RI,inplace=True)
                self.RIIS.reset_index(drop=True,inplace=True)
                print('Find %d RI stander' % (len(self.RIIS['RT'])-len(np.where(np.isnan(self.RIIS['RT'])==True)[0])))
            for i in range(1,len(self.RIIS.iloc[:,1])-1):
                second_filter = np.where(np.array(self.RIIS.at[i,'Candidate_Int'])/max(self.RIIS.at[i,'Candidate_Int'])>0.4)[0]
                if len(second_filter)>1:
                    n__1 = self.RIIS.at[i-1,'RT']
                    n_1 = self.RIIS.at[i+1,'RT']
                    score=[]
                    for ii in second_filter:
                        score.append(max(0,min(self.RIIS.at[i,'Candidate_RT'][ii]-n__1,n_1-self.RIIS.at[i,'Candidate_RT'][ii]))/max(self.RIIS.at[i,'Candidate_RT'][ii]-n__1,n_1-self.RIIS.at[i,'Candidate_RT'][ii]))
                    best = np.where(abs(np.array(score)-1)==min(abs(np.array(score)-1)))[0][0]
                    self.RIIS.at[i,'RT']=self.RIIS.at[i,'Candidate_RT'][best]
            for i in range(len(self.RIIS.iloc[:,0])):
                if i == 0:
                    if self.RIIS['RT'][i] > self.RIIS['RT'][i+1]:
                        for ii in range(self.RIIS['Candidate_RT']):
                            if ii < self.RIIS['RT'][i+1]:
                                self.RIIS['RT'][i] = ii
                                break
                        if self.RIIS['RT'][i]>self.RIIS['RT'][i+1]:
                            self.RIIS['RT'][i] = np.nan
                elif i == len(self.RIIS.iloc[:,0])-1:
                    if self.RIIS['RT'][i] < self.RIIS['RT'][i-1]:
                        for ii in range(self.RIIS['Candidate_RT']):
                            if ii > self.RIIS['RT'][i+1]:
                                self.RIIS['RT'][i] = ii
                                break
                        if self.RIIS['RT'][i] < self.RIIS['RT'][i-1]:
                            self.RIIS['RT'][i] = np.nan
                else:
                    if self.RIIS['RT'][i] < self.RIIS['RT'][i-1] or self.RIIS['RT'][i]>self.RIIS['RT'][i+1]:
                        for ii in range(self.RIIS['Candidate_RT']):
                            if ii > self.RIIS['RT'][i+1] and ii < self.RIIS['RT'][i+1]:
                                self.RIIS['RT'][i] = ii
                                break
                        if self.RIIS['RT'][i] < self.RIIS['RT'][i-1] or self.RIIS['RT'][i] > self.RIIS['RT'][i+1]:
                            self.RIIS['RT'][i] = np.nan
        if len(np.where(np.isnan(np.array(self.RIIS['RT']))==True)[0])/len(self.RIIS['RT'])>0.4:
            raise ValueError('Too much missing RI stander')
        else:
            self.RIIS.drop(self.RIIS['RT'][np.isnan(self.RIIS['RT'])==True],inplace=True)
            self.RIIS.reset_index(drop=True,inplace=True)
        self.Final_Peak_Detect['RI'] = 0
        self.Final_Peak_Detect['RI'] = self.Final_Peak_Detect['RI'].map(lambda x:('%.2f')%x)
        for i in range(len(self.Final_Peak_Detect)):
            n = self.RIIS[self.RIIS['RT']<self.Final_Peak_Detect.at[i,'RT']]
            n_1 = self.RIIS[self.RIIS['RT']>self.Final_Peak_Detect.at[i,'RT']]
            if len(n)>0 and len(n_1)>0:
                N = self.RIIS[self.RIIS['RT']<=self.Final_Peak_Detect.at[i,'RT']].index[-1]
                N_1 = self.RIIS[self.RIIS['RT']>self.Final_Peak_Detect.at[i,'RT']].index[0]
            elif len(n)>0 and len(n_1)==0:
                N = self.RIIS[self.RIIS['RT']<=self.Final_Peak_Detect.at[i,'RT']].index[-2]
                N_1 = self.RIIS[self.RIIS['RT']<=self.Final_Peak_Detect.at[i,'RT']].index[-1]
            elif len(n)==0 and len(n_1)>0:
                N = self.RIIS[self.RIIS['RT']>=self.Final_Peak_Detect.at[i,'RT']].index[0]
                N_1 = self.RIIS[self.RIIS['RT']>=self.Final_Peak_Detect.at[i,'RT']].index[1]
            RI = 100*(self.RIIS.loc[N,'C']+(self.Final_Peak_Detect.at[i,'RT']-self.RIIS.at[N,'RT'])/(self.RIIS.at[N_1,'RT']-self.RIIS.at[N,'RT']))
            self.Final_Peak_Detect.at[i,'RI'] = round(RI,2)
    def Calculate_MS2_RI(self):
        self.MS2_Data['Scan_RI'] = self.MS2_Data['Scan_Time'].apply(lambda x:EazyMZDataProcess.RT_to_RI(x,self.RIIS))
    def Output_Result(self,Path):
        mgf_output = ''
        csv_mgf = pd.DataFrame(columns=['row ID','row m/z','row retention time','correlation group ID','annotation network number','best ion','auto MS2 verify','identified by n=','partners','neutral M mass','Peak height'])
        for i in  range(len(self.Final_Peak_Detect)):
            MZ = self.Final_Peak_Detect['AverageMZ'][i]
            RT = np.around(self.Final_Peak_Detect['RT'][i]/60,3)
            Int = self.Final_Peak_Detect['Int'][i]
            temp_csv_mgf = pd.DataFrame([(i+1,MZ,RT,'','','','','','','',Int)],columns=['row ID','row m/z','row retention time','correlation group ID','annotation network number','best ion','auto MS2 verify','identified by n=','partners','neutral M mass','Peak height'])
            csv_mgf = pd.concat([csv_mgf,temp_csv_mgf])
            mgf_output = mgf_output+'BEGIN IONS\nFEATURE_ID='+str(i+1)+'\nPEPMASS='+str(MZ)+'\nSCANS='+str(i+1)+'\nRTINSECONDS='+str(RT)+'\nCHARGE=1+\nMSLEVEL=2\n'
            MS2_place = list(filter(lambda x:abs(self.MS2_Pre[x]-MZ)/MZ<0.000010 and self.Final_Peak_Detect['RTList'][i][0]<self.MS2_RT_List[x]<self.Final_Peak_Detect['RTList'][i][-1],range(len(self.MS2_Pre))))
            if len(MS2_place)>0:
                MS2_place = list(filter(lambda x:np.min(abs(self.MS2_RT_List[MS2_place]-RT))==abs(self.MS2_RT_List[x]-RT),MS2_place))[0]
                for ii in range(len(self.MS2_MZ_List[MS2_place])):
                    mgf_output = mgf_output+str(self.MS2_MZ_List[MS2_place][ii])+' '+str(self.MS2_Int_List[MS2_place][ii])+'\n'
            mgf_output = mgf_output+'END IONS\n\n'
        def namestr(obj, namespace):
            return [name for name in namespace if namespace[name] is obj]    
        with open(Path+'/'+namestr(self,globals())[0]+'.mgf','w')as mgfFile:
            mgfFile.write(mgf_output)
        csv_mgf.to_csv(Path +'/'+ namestr(self,globals())[0] +'.csv',index=False)
        OutputPath = Path +'/'+ namestr(self,globals())[0] +'.xlsx'
        self.Final_Peak_Detect['IntList'] = self.Final_Peak_Detect['IntList'].apply(lambda x:x.tolist())
        self.Final_Peak_Detect['RTList'] = self.Final_Peak_Detect['RTList'].apply(lambda x:x.tolist())
        self.Final_Peak_Detect['MS2_MZ'] = self.Final_Peak_Detect['MS2_MZ'].apply(lambda x:x.tolist())
        self.Final_Peak_Detect['MS2_Int'] = self.Final_Peak_Detect['MS2_Int'].apply(lambda x:x.tolist())
        self.Final_Peak_Detect.to_excel(OutputPath,index=False)
    def OutputSinglePeakList(self):
        Data = self.Final_Peak_Detect.copy()
        Data.sort_values(by='AverageMZ',ascending=True,inplace=True,ignore_index=True)
        Data.insert(0,'ID',[0]*len(Data))
        for i in range(len(Data)):
            Data.at[i,'ID'] = i+1
            Data.at[i,'AverageMZ'] = np.around(Data.at[i,'AverageMZ'],5)
            Data.at[i,'RT'] = np.around(Data.at[i,'RT'],2)
            Data.at[i,'RI'] = np.around(Data.at[i,'RI'],1)
        Data.drop('RTList',axis=1,inplace=True)
        Data.rename(columns={'AverageMZ':'m/z'},inplace=True)
        Data.rename(columns={'Int':'Intensity'},inplace=True)
        Data.rename(columns={'RT':'RT (s)'},inplace=True)
        Data.to_excel(self.OutputSinglePeakListPath,index=False)
    def load_FPD(self,Path):
        with open(Path,'rb') as f:
            self.Final_Peak_Detect = pickle.load(f)
    
    def save_FPD(self,Path):
        with open(Path,'wb') as f:
            pickle.dump(self.Final_Peak_Detect,f)
        

class DataAlignment(object):
    def __init__(self,Label_process_sub,process_bar):
        self.Label_process_sub = Label_process_sub
        self.process_bar = process_bar
        self.DataBase = pd.DataFrame(columns=['Data','Final_Peak_Detect','Data_Name','Tag'])
        self.AlignmentParam = {'MZ_Tor':0.000010,'RT_Tor':6,'A':0.75,
                               'Miss_Filter':0.8,'Threshold':3,'RI_Alignment':False,
                               'RT_min_Tor':6,'RI_Tor':0.02}
        self.RefList=pd.DataFrame(columns=['m/z','RT','MS_List','RT_List'])
        self.RefList['m/z'] = self.RefList['m/z'].map(lambda x:'%.4f'%x)
    def add_Data(self,Data,Tag='Sample'): #Gradient=0
        # 添加数据
        if Tag == 'Sample' or  Tag == 'QC':
            def namestr(obj, namespace):
                return [name for name in namespace if namespace[name] is obj]
            self.DataBase.at[len(self.DataBase),'Data']= ''
            self.DataBase.at[len(self.DataBase)-1,'Final_Peak_Detect']= Data['Final_Peak_Detect']
            self.DataBase.at[len(self.DataBase)-1,'Data_Name'] = Data['Name']
            self.DataBase.at[len(self.DataBase)-1,'Tag'] = Tag
            self.DataBase.sort_values(by='Tag',ascending=False,inplace=True)
            self.DataBase.reset_index(drop=True,inplace=True)
        elif Tag == 'Blank':
            self.DataBase.at[len(self.DataBase),'Data']= Data['Data']
            self.DataBase.at[len(self.DataBase)-1,'Final_Peak_Detect'] = ''
            self.DataBase.at[len(self.DataBase)-1,'Data_Name'] = Data['Data'].DataName
            self.DataBase.at[len(self.DataBase)-1,'Tag'] = Tag
            self.DataBase.sort_values(by='Tag',ascending=False,inplace=True)
            self.DataBase.reset_index(drop=True,inplace=True)
        else:
            print('Data Tag should be Sample or Blank or QC')
    def show_Data(self):
        print(self.DataBase.iloc[:,[1,2]])
    def get_param(self):
        print(self.AlignmentParam)
    def set_param(self,Name,Value):
        self.AlignmentParam[Name] = Value
    def RenewRefList(self):
        ''' 已经修改了，与FPD直接关联 '''
        # 最优权值匹配模型
        self.RefList=pd.DataFrame(columns=['m/z','Int','RT','RI','MS_List','RT_List','RI_List','MS2_Int','MS2_MZ'])
        self.RefList['m/z'] = self.RefList['m/z'].map(lambda x:'%.4f'%x)
        # 进行Alignment，获得列表
        def MergeMS2(Int_a,MZ_a,Int_b,MZ_b,MZ_Tor = self.AlignmentParam['MZ_Tor']):
            for i in range(len(Int_b)):
                same_place = np.where(abs(MZ_a-MZ_b[i])/MZ_b[i]<MZ_Tor)[0]
                if len(same_place)>0:
                    if Int_b[i]>Int_a[same_place[0]]:
                        Int_a[same_place[0]] = Int_b[i]
                        MZ_a[same_place[0]] = MZ_b[i]
                if len(same_place)==0:
                    Int_a = np.append(Int_a,Int_b[i])
                    MZ_a = np.append(MZ_a,MZ_b[i])
            return Int_a,MZ_a
        def GrubbsTest(TestList):
            ResultList = list(map(lambda x:abs(x-np.mean(TestList))/np.std(TestList),TestList))
            if len(TestList)==3:
                if max(ResultList) > 1.15:                    
                    return max(ResultList)
                else:
                    return 0
            elif len(TestList)==4:
                if max(ResultList) > 1.46:                    
                    return max(ResultList)
                else:
                    return 0
            else:
                return 0
        for i in range(len(self.DataBase)):
            self.RefList[self.DataBase.at[i,'Data_Name']] = 0
        for i in range(len(self.DataBase)):
            self.Label_process_sub.setText('Align '+str(int((i+1)/len(self.DataBase))*100)+ ' %')
            if self.DataBase.at[i,'Tag'] != 'Blank':
                self.RefList.sort_values('m/z',ascending=(True),inplace=True)
                self.RefList.reset_index(drop=True,inplace=True)
                temp_Data = self.DataBase.at[i,'Final_Peak_Detect']
                temp_Data.sort_values('AverageMZ',ascending=(True),inplace=True,ignore_index=True)
                init_RefMZ = pd.DataFrame(columns=['m/z','MZList','RowIndex'])
                init_SampleMZ = pd.DataFrame(columns=['m/z','MZList','RowIndex'])
                for ii in range(len(self.RefList)):
                    if len(init_RefMZ)==0:
                        init_RefMZ.loc[len(init_RefMZ)] = [self.RefList.at[ii,'m/z'],[self.RefList.at[ii,'m/z']],[ii]]
                    else:
                        if abs(init_RefMZ.at[len(init_RefMZ)-1,'m/z']-self.RefList.at[ii,'m/z'])/self.RefList.at[ii,'m/z']<self.AlignmentParam['MZ_Tor']:
                            init_RefMZ.at[len(init_RefMZ)-1,'MZList'].append(self.RefList.at[ii,'m/z'])
                            init_RefMZ.at[len(init_RefMZ)-1,'RowIndex'].append(ii)
                            init_RefMZ.at[len(init_RefMZ)-1,'m/z']=np.mean(init_RefMZ.at[len(init_RefMZ)-1,'MZList'])
                        else:
                            init_RefMZ.loc[len(init_RefMZ)] = [self.RefList.at[ii,'m/z'],[self.RefList.at[ii,'m/z']],[ii]]
                for ii in range(len(temp_Data)):
                    if len(init_SampleMZ)==0:
                        init_SampleMZ.loc[len(init_SampleMZ)] = [temp_Data.at[ii,'AverageMZ'],[temp_Data.at[ii,'AverageMZ']],[ii]]
                    else:
                        if abs(init_SampleMZ.at[len(init_SampleMZ)-1,'m/z']-temp_Data.at[ii,'AverageMZ'])/temp_Data.at[ii,'AverageMZ']<self.AlignmentParam['MZ_Tor']:
                            init_SampleMZ.at[len(init_SampleMZ)-1,'MZList'].append(temp_Data.at[ii,'AverageMZ'])
                            init_SampleMZ.at[len(init_SampleMZ)-1,'RowIndex'].append(ii)
                            init_SampleMZ.at[len(init_SampleMZ)-1,'m/z']=np.mean(init_SampleMZ.at[len(init_SampleMZ)-1,'MZList'])
                        else:
                            init_SampleMZ.loc[len(init_SampleMZ)] = [temp_Data.at[ii,'AverageMZ'],[temp_Data.at[ii,'AverageMZ']],[ii]]
                add_MZ = []
                add_Int = []
                add_RT = []
                add_RI = []
                add_MS_List = []
                add_RT_List = [] 
                add_RI_List = []
                if self.AlignmentParam['RI_Alignment'] == False:
                    Sample_RT = temp_Data.loc[:,'RT']
                    init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True)
                    self.process_bar.setMaximum(len(init_SampleMZ))
                    self.process_bar.setValue(0)
                    for ii in range(len(init_SampleMZ)):
                        self.process_bar.setValue(self.process_bar.value()+1)
                        QApplication.processEvents()
                        Sample_MZ = init_SampleMZ.at[ii,'m/z']
                        #match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
                        match_Index = range(bisect.bisect_left(init_RefMZ['m/z'],Sample_MZ*(1-self.AlignmentParam['MZ_Tor'])),bisect.bisect_right(init_RefMZ['m/z'],Sample_MZ*(1+self.AlignmentParam['MZ_Tor'])))
                        if len(match_Index)>1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                            for i_mI in range(1,len(match_Index)):
                                Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                        elif len(match_Index)==1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                        else:
                            for i_add in init_SampleMZ.at[ii,'RowIndex']:
                                add_MZ.append(temp_Data.at[i_add,'AverageMZ'])
                                add_Int.append(temp_Data.at[i_add,'Int'])
                                add_RT.append(Sample_RT[i_add])
                                add_RI.append(temp_Data.at[i_add,'RI'])
                                add_MS_List.append([temp_Data.at[i_add,'AverageMZ']])
                                add_RT_List.append([Sample_RT[i_add]])
                                add_RI_List.append([temp_Data.at[i_add,'RI']])
                            continue
                        Score_matrix= np.zeros([len(Ref_Index),len(init_SampleMZ.at[ii,'RowIndex'])])
                        for i_Ref in range(len(Ref_Index)):
                            for i_Sample in range(len(init_SampleMZ.at[ii,'RowIndex'])):
                                i_Ref_MZ = self.RefList.at[Ref_Index[i_Ref],'m/z']
                                i_Sample_MZ = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sample],'AverageMZ']
                                i_Ref_Time = self.RefList.at[Ref_Index[i_Ref],'RT']
                                i_Sample_Time = Sample_RT[init_SampleMZ.at[ii,'RowIndex'][i_Sample]]
                                if abs(i_Ref_MZ-i_Sample_MZ)/i_Ref_MZ<self.AlignmentParam['MZ_Tor'] and abs(i_Ref_Time-i_Sample_Time)<self.AlignmentParam['RT_Tor']:
                                    Score_matrix[i_Ref,i_Sample] = self.AlignmentParam['A']*np.exp(-0.5*((i_Sample_Time-i_Ref_Time)/(self.AlignmentParam['RT_Tor']))**2)+(1-self.AlignmentParam['A'])*np.exp(-0.5*((i_Sample_MZ-i_Ref_MZ)/(i_Sample_MZ*self.AlignmentParam['MZ_Tor']))**2)
                        # 总体最优值匹配
                        Sm_row,Sm_col = linear_sum_assignment(Score_matrix,True)
                        for i_Sm_row,i_Sm_col in zip(Sm_row,Sm_col):
                            if Score_matrix[i_Sm_row,i_Sm_col] != 0:
                                self.RefList.at[Ref_Index[i_Sm_row],self.DataBase.at[i,'Data_Name']] = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'Int']
                                self.RefList.at[Ref_Index[i_Sm_row],'MS_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                                self.RefList.at[Ref_Index[i_Sm_row],'RT_List'].append(Sample_RT[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                                self.RefList.at[Ref_Index[i_Sm_row],'RI_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'RI'])
                                self.RefList.at[Ref_Index[i_Sm_row],'RT'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'RT_List'])
                                self.RefList.at[Ref_Index[i_Sm_row],'RI'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'RI_List'])
                                self.RefList.at[Ref_Index[i_Sm_row],'m/z'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'MS_List'])
                            else:
                                add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                                add_Int.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'Int'])
                                add_RT.append(Sample_RT[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                                add_RI.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'RI'])
                                add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ']])
                                add_RT_List.append([Sample_RT[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]]])
                                add_RI_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'RI']])
                        miss_col = list(filter(lambda x:x not in Sm_col,range(len(init_SampleMZ.at[ii,'RowIndex']))))
                        for i_miss in miss_col:
                            add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ'])
                            add_Int.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'Int'])
                            add_RT.append(Sample_RT[init_SampleMZ.at[ii,'RowIndex'][i_miss]])
                            add_RI.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'RI'])
                            add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ']])
                            add_RT_List.append([Sample_RT[init_SampleMZ.at[ii,'RowIndex'][i_miss]]])
                            add_RI_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'RI']])
                else:
                    Sample_RI = temp_Data.loc[:,'RI']
                    init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True)
                    self.process_bar.setMaximum(len(init_SampleMZ))
                    self.process_bar.setValue(0)
                    #bar = Bar('Alignment '+str(i+1)+' / '+str(len(self.DataBase)), max=len(init_SampleMZ))
                    for ii in range(len(init_SampleMZ)):
                        #bar.next()
                        self.process_bar.setValue(self.process_bar.value()+1)
                        QApplication.processEvents()
                        Sample_MZ = init_SampleMZ.at[ii,'m/z']
                        match_Index = range(bisect.bisect_left(init_RefMZ['m/z'],Sample_MZ*(1-self.AlignmentParam['MZ_Tor'])),bisect.bisect_right(init_RefMZ['m/z'],Sample_MZ*(1+self.AlignmentParam['MZ_Tor'])))
                        #match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
                        if len(match_Index)>1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                            for i_mI in range(1,len(match_Index)):
                                Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                        elif len(match_Index)==1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                        else:
                            for i_add in init_SampleMZ.at[ii,'RowIndex']:
                                add_MZ.append(temp_Data.at[i_add,'AverageMZ'])
                                add_Int.append(temp_Data.at[i_add,'Int'])
                                add_RT.append(temp_Data.at[i_add,'RT'])
                                add_RI.append(Sample_RI[i_add])
                                add_MS_List.append([temp_Data.at[i_add,'AverageMZ']])
                                add_RT_List.append([temp_Data.at[i_add,'RT']])
                                add_RI_List.append([Sample_RI[i_add]])
                            continue
                        Score_matrix= np.zeros([len(Ref_Index),len(init_SampleMZ.at[ii,'RowIndex'])])
                        for i_Ref in range(len(Ref_Index)):
                            for i_Sample in range(len(init_SampleMZ.at[ii,'RowIndex'])):
                                i_Ref_MZ = self.RefList.at[Ref_Index[i_Ref],'m/z']
                                i_Sample_MZ = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sample],'AverageMZ']
                                i_Ref_Time = self.RefList.at[Ref_Index[i_Ref],'RI']
                                i_Sample_Time = Sample_RI[init_SampleMZ.at[ii,'RowIndex'][i_Sample]]
                                if abs(i_Ref_MZ-i_Sample_MZ)/i_Ref_MZ<self.AlignmentParam['MZ_Tor'] and abs(i_Ref_Time-i_Sample_Time)/i_Ref_Time<self.AlignmentParam['RI_Tor']:
                                    Score_matrix[i_Ref,i_Sample] = self.AlignmentParam['A']*np.exp(-0.5*((i_Sample_Time-i_Ref_Time)/(i_Ref_Time*self.AlignmentParam['RI_Tor']))**2)+(1-self.AlignmentParam['A'])*np.exp(-0.5*((i_Sample_MZ-i_Ref_MZ)/(i_Sample_MZ*self.AlignmentParam['MZ_Tor']))**2)
                        # 总体最优值匹配
                        Sm_row,Sm_col = linear_sum_assignment(Score_matrix,True)
                        for i_Sm_row,i_Sm_col in zip(Sm_row,Sm_col):
                            if Score_matrix[i_Sm_row,i_Sm_col] != 0:
                                self.RefList.at[Ref_Index[i_Sm_row],self.DataBase.at[i,'Data_Name']] = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'Int']
                                self.RefList.at[Ref_Index[i_Sm_row],'MS_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                                self.RefList.at[Ref_Index[i_Sm_row],'RT_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'RT'])
                                self.RefList.at[Ref_Index[i_Sm_row],'RI_List'].append(Sample_RI[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                                self.RefList.at[Ref_Index[i_Sm_row],'RI'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'RI_List'])
                                self.RefList.at[Ref_Index[i_Sm_row],'RT'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'RT_List'])
                                self.RefList.at[Ref_Index[i_Sm_row],'m/z'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'MS_List'])
                            else:
                                add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                                add_Int.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'Int'])
                                add_RI.append(Sample_RI[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                                add_RT.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'RT'])
                                add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ']])
                                add_RI_List.append([Sample_RI[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]]])
                                add_RT_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'RT']])
                        miss_col = list(filter(lambda x:x not in Sm_col,range(len(init_SampleMZ.at[ii,'RowIndex']))))
                        for i_miss in miss_col:
                            add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ'])
                            add_Int.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'Int'])
                            add_RT.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'RT'])
                            add_RI.append(Sample_RI[init_SampleMZ.at[ii,'RowIndex'][i_miss]])
                            add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ']])
                            add_RT_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'RT']])
                            add_RI_List.append([Sample_RI[init_SampleMZ.at[ii,'RowIndex'][i_miss]]])
                    #bar.finish()
                
                add_RefList = pd.DataFrame({'m/z': add_MZ,'RT': add_RT,'RI':add_RI,self.DataBase.at[i, 'Data_Name']: add_Int,
                                            'MS_List':add_MS_List,'RT_List':add_RT_List,'RI_List':add_RI_List})
                self.RefList = pd.concat([self.RefList, add_RefList])
                self.RefList.reset_index(drop=True, inplace=True)
                self.RefList = self.RefList.fillna(0)
        self.Label_process_sub.setText('Blank Filter')
        QApplication.processEvents()
        Sample_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Sample'])
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        self.RefList['max_SampleInt'] = self.RefList.apply(lambda x:max(x[Sample_Name]),axis=1)
        self.RefList['Int'] = self.RefList.apply(lambda x:np.mean(x[Sample_Name]),axis=1)
        if len(Blank_Name) > 0:
            for i in Blank_Name:
                temp_Blank_Data = self.DataBase[self.DataBase['Data_Name']==i].copy()
                temp_Blank_Data.reset_index(drop=True,inplace=True)
                temp_Blank_Data = temp_Blank_Data['Data'][0]
                self.process_bar.setMaximum(len(self.RefList))
                self.process_bar.setValue(0)
                for ii in range(len(self.RefList)):
                    self.process_bar.setValue(self.process_bar.value()+1)
                    QApplication.processEvents()
                    [Blank_RT_List,Blank_Int_List] = temp_Blank_Data.ExtractDataPoint(self.RefList.at[ii,'m/z'],self.RefList.at[ii,'RT']+3,self.RefList.at[ii,'RT']-3,smooth_index=0)
                    self.RefList.loc[ii,i] = max(1,max(Blank_Int_List))
            self.RefList['mean_BlankInt'] = self.RefList.apply(lambda x:np.mean(x[Blank_Name]) if max(x[Blank_Name])>0 else 1,axis=1)
        self.Label_process_sub.setText('Align finished')
    
    def assign_MS2_KDtree(self,MergeMS2,temp_RI_Tor,temp_MZ_Tor,V_List,RI_swich=False):
        self.RefList['MS2_MZ'] = self.RefList['m/z'].apply(lambda x:[])
        self.RefList['MS2_Int'] = self.RefList['m/z'].apply(lambda x:[])
        self.RefList['MS2_MZ'] = self.RefList['MS2_MZ'].astype('object')
        self.RefList['MS2_Int'] = self.RefList['MS2_Int'].astype('object')
        if RI_swich == True:
            data_points = np.column_stack((MergeMS2['Scan_RI'],MergeMS2['Pre_MZ']))
            kdtree = KDTree(data_points)
            self.Label_process_sub.setText('Assign MS2')
            self.process_bar.setMaximum(len(self.RefList))
            self.process_bar.setValue(0)
            for i,RI,MZ in zip(range(len(self.RefList)),self.RefList['RI'],self.RefList['m/z']):
                self.process_bar.setValue(self.process_bar.value()+1)
                QApplication.processEvents()
                radius = RI*temp_RI_Tor
                Match_List = kdtree.query_ball_point([RI,MZ], r=radius)
                Match_List = list(filter(lambda x:abs(MergeMS2.at[x,'Pre_MZ']-MZ)/MZ<temp_MZ_Tor,Match_List))
                if len(Match_List)>0:
                    Match_List = list(filter(lambda x:abs(RI-MergeMS2.loc[x,'Scan_RI'])==min(abs(np.array(MergeMS2.loc[Match_List,'Scan_RI'])-RI)),Match_List))[0]
                    CE_Index = list(filter(lambda x:V_List[x]==MergeMS2.at[Match_List,'CE'],range(len(V_List))))[0]
                    if Match_List-CE_Index > 0 and Match_List+len(V_List)-CE_Index < len(MergeMS2):
                        Valified_Index_List = list(range(Match_List-CE_Index,Match_List+len(V_List)-CE_Index))
                    elif Match_List-CE_Index < 0:
                        Valified_Index_List = list(range(0,Match_List+len(V_List)-CE_Index))
                    elif Match_List+len(V_List)-CE_Index > len(MergeMS2):
                        Valified_Index_List = list(range(Match_List-CE_Index,len(MergeMS2)))
                    MZ_List = []
                    Int_List = []
                    for ii in Valified_Index_List:
                        if abs(MergeMS2.at[ii,'Pre_MZ'] - MZ) / MZ<temp_MZ_Tor*2 and len(MergeMS2.at[ii,'Int_List'])>0:
                            if len(MZ_List) == 0:
                                MZ_List = list(MergeMS2.at[ii,'MZ_List'][MergeMS2.at[ii,'Int_List']>max(MergeMS2.at[ii,'Int_List'])*0.05])
                                Int_List = list(MergeMS2.at[ii,'Int_List'][MergeMS2.at[ii,'Int_List']>max(MergeMS2.at[ii,'Int_List'])*0.05])
                            else:
                                #min_limit = max(MergeMS2.at[ii,'Int_List'])*0.01
                                for iii in range(len(MergeMS2.at[ii,'MZ_List'])):
                                    #if MergeMS2.at[ii,'Int_List'][iii] >= min_limit:
                                        MS2_Merge_Index = list(filter(lambda x:abs(MZ_List[x]-MergeMS2.at[ii,'MZ_List'][iii])/MZ_List[x]<temp_MZ_Tor,range(len(MZ_List))))
                                        if len(MS2_Merge_Index) == 0:
                                            MZ_List.append(MergeMS2.at[ii,'MZ_List'][iii])
                                            Int_List.append(MergeMS2.at[ii,'Int_List'][iii])
                                        elif Int_List[MS2_Merge_Index[0]] < MergeMS2.at[ii,'Int_List'][iii]:
                                            MZ_List[MS2_Merge_Index[0]] = MergeMS2.at[ii,'MZ_List'][iii]
                                            Int_List[MS2_Merge_Index[0]] = MergeMS2.at[ii,'Int_List'][iii]
                    temp_range = pd.DataFrame({'MZ_List':MZ_List,'Int_List':Int_List})
                    temp_range.sort_values(by='MZ_List',inplace=True)
                    temp_range.reset_index(drop=True,inplace=True)
                    self.RefList.at[i,'MS2_MZ'] = list(temp_range['MZ_List'])
                    self.RefList.at[i,'MS2_Int'] = list(temp_range['Int_List'])
        elif RI_swich == False:
            data_points = np.column_stack((MergeMS2['Scan_Time'],MergeMS2['Pre_MZ']))
            kdtree = KDTree(data_points)
            radius = temp_RI_Tor
            self.Label_process_sub.setText('Assign MS2')
            self.process_bar.setMaximum(len(self.RefList))
            self.process_bar.setValue(0)
            for i,RT,MZ in zip(range(len(self.RefList)),self.RefList['RT'],self.RefList['m/z']):
                self.process_bar.setValue(self.process_bar.value()+1)
                QApplication.processEvents()
                Match_List = kdtree.query_ball_point([RT,MZ], r=radius)
                Match_List = list(filter(lambda x:abs(MergeMS2.at[x,'Pre_MZ']-MZ)/MZ<temp_MZ_Tor,Match_List))
                if len(Match_List)>0:
                    Match_List = list(filter(lambda x:abs(RT-MergeMS2.loc[x,'Scan_Time'])==min(abs(np.array(MergeMS2.loc[Match_List,'Scan_Time'])-RT)),Match_List))[0]
                    CE_Index = list(filter(lambda x:V_List[x]==MergeMS2.at[Match_List,'CE'],range(len(V_List))))[0]
                    if Match_List-CE_Index > 0 and Match_List+len(V_List)-CE_Index < len(MergeMS2):
                        Valified_Index_List = list(range(Match_List-CE_Index,Match_List+len(V_List)-CE_Index))
                    elif Match_List-CE_Index < 0:
                        Valified_Index_List = list(range(0,Match_List+len(V_List)-CE_Index))
                    elif Match_List+len(V_List)-CE_Index > len(MergeMS2):
                        Valified_Index_List = list(range(Match_List-CE_Index,len(MergeMS2)))
                    MZ_List = []
                    Int_List = []
                    for ii in Valified_Index_List:
                        if abs(MergeMS2.at[ii,'Pre_MZ'] - MZ) / MZ<temp_MZ_Tor*2 and len(MergeMS2.at[ii,'Int_List'])>0:
                            if len(MZ_List) == 0:
                                MZ_List = list(MergeMS2.at[ii,'MZ_List'][MergeMS2.at[ii,'Int_List']>max(MergeMS2.at[ii,'Int_List'])*0.05])
                                Int_List = list(MergeMS2.at[ii,'Int_List'][MergeMS2.at[ii,'Int_List']>max(MergeMS2.at[ii,'Int_List'])*0.05])
                            else:
                                #min_limit = max(MergeMS2.at[ii,'Int_List'])*0.01
                                for iii in range(len(MergeMS2.at[ii,'MZ_List'])):
                                    #if MergeMS2.at[ii,'Int_List'][iii] >= min_limit:
                                        MS2_Merge_Index = list(filter(lambda x:abs(MZ_List[x]-MergeMS2.at[ii,'MZ_List'][iii])/MZ_List[x]<self.get_param('MS1_Tor'),range(len(MZ_List))))
                                        if len(MS2_Merge_Index) == 0:
                                            MZ_List.append(MergeMS2.at[ii,'MZ_List'][iii])
                                            Int_List.append(MergeMS2.at[ii,'Int_List'][iii])
                                        elif Int_List[MS2_Merge_Index[0]] < MergeMS2.at[ii,'Int_List'][iii]:
                                            MZ_List[MS2_Merge_Index[0]] = MergeMS2.at[ii,'MZ_List'][iii]
                                            Int_List[MS2_Merge_Index[0]] = MergeMS2.at[ii,'Int_List'][iii]
                    temp_range = pd.DataFrame({'MZ_List':MZ_List,'Int_List':Int_List})
                    temp_range.sort_values(by='MZ_List',inplace=True)
                    temp_range.reset_index(drop=True,inplace=True)
                    self.RefList.at[i,'MS2_MZ'] = list(temp_range['MZ_List'])
                    self.RefList.at[i,'MS2_Int'] = list(temp_range['Int_List'])

    def pool_MS2_assignment(MZ,RI,MergeMS2,temp_RI_Tor,temp_MZ_Tor,V_List):
        #Valified_List = list(filter(lambda x:abs(MZ-MergeMS2.at[x,'Pre_MZ'])/MZ<temp_MZ_Tor and RI*(1-temp_RI_Tor)<MergeMS2.at[x,'Scan_RI']<RI*((1+temp_RI_Tor)),range(len(MergeMS2))))
        Valified_List = list(filter(lambda x:abs(MZ-MergeMS2.at[x,'Pre_MZ'])/MZ<temp_MZ_Tor,range(len(MergeMS2))))
        if len(Valified_List) >0 :
            Valified_List = list(filter(lambda x:RI*(1-temp_RI_Tor)<MergeMS2.at[x,'Scan_RI']<RI*((1+temp_RI_Tor)),Valified_List))
            if len(Valified_List) >0 :
                Valified_Index = Valified_List[np.where(abs(np.array(MergeMS2.loc[Valified_List,'Scan_RI'])-RI)==min(abs(np.array(MergeMS2.loc[Valified_List,'Scan_RI'])-RI)))[0][0]]
                CE_Index = list(filter(lambda x:V_List[x]==MergeMS2.at[Valified_Index+3,'CE'],range(len(V_List))))[0]
                if Valified_Index-CE_Index > 0 and Valified_Index+len(V_List)-CE_Index < len(MergeMS2):
                    Valified_Index_List = list(range(Valified_Index-CE_Index,Valified_Index+len(V_List)-CE_Index))
                elif Valified_Index-CE_Index < 0:
                    Valified_Index_List = list(range(0,Valified_Index+len(V_List)-CE_Index))
                elif Valified_Index+len(V_List)-CE_Index > len(MergeMS2):
                    Valified_Index_List = list(range(Valified_Index-CE_Index,len(MergeMS2)))
                MZ_List = []
                Int_List = []
                for ii in Valified_Index_List:
                    if abs(MergeMS2.at[ii,'Pre_MZ'] - MZ) / MZ<temp_MZ_Tor:
                        if len(MZ_List) == 0:
                            MZ_List = list(MergeMS2.at[ii,'MZ_List'])
                            Int_List = list(MergeMS2.at[ii,'Int_List'])
                        else:
                            Threshold = max(MergeMS2.at[ii,'Int_List'])*0.05
                            for iii in range(len(MergeMS2.at[ii,'MZ_List'])):
                                if MergeMS2.at[ii,'Int_List'][iii] > Threshold:
                                    MS2_Merge_Index = list(filter(lambda x:abs(MZ_List[x]-MergeMS2.at[ii,'MZ_List'][iii])/MZ_List[x]<temp_MZ_Tor,range(len(MZ_List))))
                                    if len(MS2_Merge_Index) == 0:
                                        MZ_List.append(MergeMS2.at[ii,'MZ_List'][iii])
                                        Int_List.append(MergeMS2.at[ii,'Int_List'][iii])
                                    elif Int_List[MS2_Merge_Index[0]] < MergeMS2.at[ii,'Int_List'][iii]:
                                        MZ_List[MS2_Merge_Index[0]] = MergeMS2.at[ii,'MZ_List'][iii]
                                        Int_List[MS2_Merge_Index[0]] = MergeMS2.at[ii,'Int_List'][iii]
                temp_range = pd.DataFrame({'MZ_List':MZ_List,'Int_List':Int_List})
                temp_range.sort_values(by='MZ_List',inplace=True,ignore_index=True)
                #temp_range.reset_index(drop=True,inplace=True)
                MZ_List = list(temp_range['MZ_List'])
                Int_List = list(temp_range['Int_List'])
                return [MZ_List,Int_List]
            else:
                return 0
        else:
            return 0
            
    def BlankFilter(self):
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        if len(Blank_Name) > 0:
            self.RefList.drop(self.RefList[self.RefList['max_SampleInt']/self.RefList['mean_BlankInt']<self.AlignmentParam['Threshold']].index,inplace=True)
            self.RefList.reset_index(drop=True,inplace=True)
    def Filter_MissingValue(self,WhetherDel=True):
        bar = Bar('Missing Values Filter',max = len(self.RefList))
        Sample_Number = list(filter(lambda x:self.DataBase.at[x,'Tag']=='Sample',range(len(self.DataBase))))
        Data_Name = list(self.DataBase.loc[Sample_Number,'Data_Name'])
        self.RefList['Count']=0
        self.RefList['Count'] = self.RefList['Count'].map(lambda x:'%.2f'%x)
        del_row = []
        for i in range(len(self.RefList)):
            bar.next()
            temp_row = list(filter(lambda x:self.RefList.at[i,x]!=0,Data_Name))
            self.RefList.at[i,'Count'] = len(temp_row)/len(Data_Name)
            if len(temp_row)/len(Data_Name) < self.AlignmentParam['Miss_Filter']:
                del_row.append(i)
        bar.finish()
        if WhetherDel==True:
            self.RefList.drop(del_row,inplace=True)
            self.RefList.reset_index(drop=True,inplace=True)
    
    def Output_Result(self,OutputPath):
        self.RefList.to_excel(OutputPath,index=False)



if __name__ == '__main__':
    mp.freeze_support()
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    app.setWindowIcon(QIcon('./MultipleGradientProcessor.ico'))
    app.setQuitOnLastWindowClosed(True)
    #main = FirstMainWindow()
    main_HP = HeuristicProcessorUI()
    main_HP.show()
    sys.exit(app.exec_())
        