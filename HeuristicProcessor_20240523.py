#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 20:24:43 2023

@author: lacter
"""
import pandas as pd
import time
from progress.bar import Bar
from matplotlib import pyplot as plt
import pyopenms
import math
import numpy as np
import bisect
import re
import os
from scipy.optimize import linear_sum_assignment
import warnings
'''---------'''
import sys
from PyQt5.QtWidgets import * 
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import os
from PyQt5 import QtWidgets,QtCore
import pathos
import multiprocessing as mp
import pickle

class HeuristicProcessorUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
        # 设置标题
        self.resize(600,700)
        self.centerWindow()
        
        # 设置参数字典
        self.Data_params={}
        self.Calculate_params={'MS1_Tor':0.000010,'smooth':5,'min_Int':9000,'Points':17}
    
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
        self.Lable_RT_Tor = QLabel('RT Tolerance(%)')
        self.Lable_Int_min = QLabel('Min intensity')
        self.Lable_Point = QLabel('Min scan points')
        self.Lable_SN = QLabel('Signal/Noise')
        self.Lable_SB = QLabel('Signal/Blank')
        self.LineEdit_MS_Tor = QLineEdit('10')
        self.LineEdit_RT_Tor = QLineEdit('2')
        self.LineEdit_Int_min = QLineEdit('9000')
        self.LineEdit_Point = QLineEdit('17')
        self.LineEdit_SN = QLineEdit('10')
        self.LineEdit_SB = QLineEdit('5')
        self.LineEdit_MS_Tor.setToolTip('Under this threshold will be recognized as same m/z')
        self.LineEdit_RT_Tor.setToolTip('Under this threshold will be recognized as same RT')
        self.LineEdit_Int_min.setToolTip('Under this threshold will be deleted')
        self.LineEdit_Point.setToolTip('Scan points under this threshold will be deleted')
        self.LineEdit_SN.setToolTip('Under this threshold will be deleted')
        self.LineEdit_SB.setToolTip('Under this threshold will be deleted')
        Params_setting_Layout.addWidget(self.Label_Params_Select,0,0)
        #Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.LineEdit_MS_Tor,1,1)
        Params_setting_Layout.addWidget(self.Lable_RT_Tor,1,2)
        Params_setting_Layout.addWidget(self.LineEdit_RT_Tor,1,3)
        Params_setting_Layout.addWidget(self.Lable_Int_min,1,4)
        Params_setting_Layout.addWidget(self.LineEdit_Int_min,1,5)
        Params_setting_Layout.addWidget(self.Lable_Point,2,0)
        Params_setting_Layout.addWidget(self.LineEdit_Point,2,1)
        Params_setting_Layout.addWidget(self.Lable_SN,2,2)
        Params_setting_Layout.addWidget(self.LineEdit_SN,2,3)
        Params_setting_Layout.addWidget(self.Lable_SB,2,4)
        Params_setting_Layout.addWidget(self.LineEdit_SB,2,5)
        Params_setting_Widget.setLayout(Params_setting_Layout)
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
        FileName,FileType = QFileDialog.getOpenFileNames(self,"选取文件",os.getcwd(),"mzML Files(*.mzML)")
        for i in FileName:
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
                name_begin = FileName.rfind('\\')
            name_end = len(i)
            self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
            self.TextBrowser_SampleSelect.setItem(self.TextBrowser_SampleSelect.rowCount()-1,0,QTableWidgetItem(i[name_begin+1:name_end-5]))
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,Combox_Type)
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,2,LineEdit_Group)
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,3,LineEdit_Index)
            self.Data_params[i[name_begin+1:name_end-5]] = {'Name':i[name_begin+1:name_end-5],'Path':i,'Type':'Sample','Group':1,'Index':str(self.TextBrowser_SampleSelect.rowCount()),'Combox':Combox_Type,'LineEdit_Index':LineEdit_Index,'LineEdit_Group':LineEdit_Group}

    def SelectRIISFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"选取文件",os.getcwd(),"Excel Files(*.xlsx)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_RIISSelect.setText(FileName[name_begin+1:name_end])
        self.RIISPath = FileName
        
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
    
    def Run(self):
        '''
        for i in self.Data_params.keys():
            pass
        '''
        self.SampleData = {}
        self.BlankData = {}
        self.SampleClass = {}
        self.MS2_Data = {}
        for i in self.Data_params.keys():
            if self.Data_params[i]['Type'] == 'Sample':
                self.SampleData[i] = self.Data_params[i]
                self.SampleClass[i] = self.Data_params[i]['Group']
            elif self.Data_params[i]['Type'] == 'MS2':
                self.MS2_Data[i] = self.Data_params[i]
            elif self.Data_params[i]['Type'] == 'Blank':
                self.BlankData[i] = self.Data_params[i]
        for i in self.SampleData.keys():
            if sys.platform == 'darwin':
                name_begin = self.SampleData[i]['Path'].rfind('/')
                self.filepath_title = self.SampleData[i]['Path'][0:name_begin]+'/'
            elif sys.platform == 'win':
                name_begin = self.SampleData[i]['Path'].rfind('\\')
                self.filepath_title = self.SampleData[i]['Path'][0:name_begin]+'\\'
            break
        ClassList = []
        for i in self.SampleClass.keys():
            ClassList.append(self.SampleClass[i])
        ClassList = list(set(ClassList))
        self.Align_path = self.filepath_title+'Alignment-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx'
        self.Label_process_sub.setText('Load Data')
        self.process_bar.setMaximum(len(self.SampleData.keys()))
        self.process_bar.setValue(0)
        pool = mp.Pool(os.cpu_count()-2)
        pool_result = {}
        for i in self.SampleData.keys():
            print('Load Data '+str(i))
            pool_result[i] = pool.apply_async(EazyMZDataProcess,args=(self.SampleData[i]['Path'],))
        pool.close()
        pool.join()
        for i in self.SampleData.keys():
            self.SampleData[i]['Data'] = pool_result[i].get()
        for i in self.SampleData.keys():
            self.process_bar.setValue(self.process_bar.value()+1)
            QApplication.processEvents()
            self.SampleData[i]['Data'].set_param('MS1_Tor',float(self.LineEdit_MS_Tor.text())/1000000)
            self.SampleData[i]['Data'].set_param('RT_Tor',float(self.LineEdit_RT_Tor.text())/100)
            self.SampleData[i]['Data'].set_param('min_Int',int(self.LineEdit_Int_min.text()))
            self.SampleData[i]['Data'].set_param('Points',int(self.LineEdit_Point.text()))
            self.SampleData[i]['Data'].set_RIIS(self.RIISPath)
        temp_key = list(self.SampleData.keys())[0]
        self.Polarity = self.SampleData[temp_key]['Data'].get_param('Polarity')
        self.Label_process_sub.setText('Load Blank Data')
        self.process_bar.setMaximum(len(self.BlankData.keys()))
        self.process_bar.setValue(0)
        for i in self.BlankData.keys():
            self.process_bar.setValue(self.process_bar.value()+1)
            QApplication.processEvents()
            self.BlankData[i]['Data'] = EazyMZDataProcess(self.BlankData[i]['Path'])
            self.BlankData[i]['Data'].set_RIIS(self.RIISPath)
            self.BlankData[i]['Data'].set_param('MS1_Tor',float(self.LineEdit_MS_Tor.text())/1000000)
            self.BlankData[i]['Data'].set_param('RT_Tor',float(self.LineEdit_Point.text())/100)
        self.Label_process_sub.setText('Load MS2 Data')
        self.process_bar.setMaximum(len(self.MS2_Data.keys()))
        self.process_bar.setValue(0)
        for i in self.MS2_Data.keys():
            self.process_bar.setValue(self.process_bar.value()+1)
            QApplication.processEvents()
            self.MS2_Data[i]['Data'] = EazyMZDataProcess(self.MS2_Data[i]['Path'])
            self.MS2_Data[i]['Data'].set_RIIS(self.RIISPath)
            self.MS2_Data[i]['Data'].set_param('MS1_Tor',float(self.LineEdit_MS_Tor.text())/1000000)
            self.MS2_Data[i]['Data'].set_param('RT_Tor',float(self.LineEdit_Point.text())/100)
        self.Label_process_sub.setText('Set Heuristic List')
        MS_Tor = float(self.LineEdit_MS_Tor.text())/1000000
        with open(self.HeuristicListPath,'rb') as f:
            Alltemp_HPL = pickle.load(f)
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
                same_list = list(filter(lambda x:EazyMZDataProcess.ppm_compare(MZ,temp_Heuristic_peak_list.at[x,'AverageMZ'])<MS_Tor and 
                                        ((RI_left<=temp_Heuristic_peak_list.at[x,'RI_right'] and RI_right>=temp_Heuristic_peak_list.at[x,'RI_right']) or 
                                         (RI_left<=temp_Heuristic_peak_list.at[x,'RI_left'] and RI_right>=temp_Heuristic_peak_list.at[x,'RI_left'])),range(len(temp_Heuristic_peak_list))))
                contain_list = list(filter(lambda x:EazyMZDataProcess.ppm_compare(MZ,temp_Heuristic_peak_list.at[x,'AverageMZ'])<MS_Tor and (RI_left>=temp_Heuristic_peak_list.at[x,'RI_left'] and RI_right<=temp_Heuristic_peak_list.at[x,'RI_right']),range(len(temp_Heuristic_peak_list))))
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
        self.HeuristicList =  temp_Heuristic_peak_list.copy()
        for i in self.SampleData.keys():
            self.SampleData[i]['Data'].Heuristic_peak_list = self.HeuristicList.copy()
        self.Label_process_sub.setText('Heuristic detect')
        pool = mp.Pool(os.cpu_count()-2)
        pool_result = {}
        self.process_bar.setMaximum(len(self.SampleData.keys()))
        self.process_bar.setValue(0)
        QApplication.processEvents()
        for i in self.SampleData.keys():
            print('Load Data '+str(i))
            pool_result[i]=pool.apply_async(self.SampleData[i]['Data'].Heuristic_PeakDetect,callback=self.processbar_fresh)
        pool.close()
        pool.join()
        for i in self.SampleData.keys():
            self.SampleData[i]['Data'].Final_Peak_Detect = pool_result[i].get()[0].copy()
            self.SampleData[i]['Data'].Heuristic_peak_list = pool_result[i].get()[1].copy()
        self.Label_process.setText('Calculate RI & SN')
        for i in self.SampleData.keys():
            self.SampleData[i]['Data'].Calculate_RI()
            self.SampleData[i]['Data'].Calculate_SN(drop=True,Threshold=int(self.LineEdit_SN.text()))
        pool = mp.Pool(os.cpu_count()-2)
        for i in self.SampleData.keys():
            pool.apply_async(self.SampleData[i]['Data'].Final_Peak_Detect.to_excel,(self.SampleData[i]['Path'][0:len(self.SampleData[i]['Path'])-5]+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',))
        pool.close()
        pool.join()
        '''
        for i in self.SampleData.keys():
            with open(self.SampleData[i]['Path'][0:len(self.SampleData[i]['Path'])-5]+'.pkl','wb') as f:
                f.write(pickle.dumps(self.SampleData[i]['Data']))
        with open(self.SampleData[list(main.SampleData.keys())[0]]['Path'][0:name_begin+1]+'DataStore'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl','wb') as f:
            all_data = {}
            for i in self.SampleData.keys():
                all_data[i] = self.SampleData[i]['Path'][0:len(self.SampleData[i]['Path'])-5]+'.pkl'
            f.write(pickle.dumps(all_data))
        '''
        self.Label_process.setText('Data Alignment')
        self.Align = DataAlignment(self.Label_process,self.process_bar)
        for i in self.SampleData.keys():
            self.Align.add_Data(self.SampleData[i]['Data'],Tag='Sample')
        for i in self.BlankData.keys():
            self.Align.add_Data(self.BlankData[i]['Data'],Tag='Blank')
        self.Align.set_param('RI_Alignment',True)
        self.Align.set_param('MZ_Tor',float(self.LineEdit_MS_Tor.text())/1000000)
        self.Align.set_param('RT_Tor',float(self.LineEdit_Point.text())/100)
        self.Align.RenewRefList()
        ''' MS2 assine'''
        for i in self.MS2_Data.keys():  
            V_List_Num = list(self.MS2_Data[i]['Data'].MS2_Data['Pre_MZ'][int(len(self.MS2_Data[i]['Data'].MS2_Data['Pre_MZ'])/2):int(len(self.MS2_Data[i]['Data'].MS2_Data['Pre_MZ'])/2)+50])
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
            self.MS2_Data[i]['Data'].MS2_Data['CE'] = V_List[0]
            for ii in range(1,len(self.MS2_Data[i]['Data'].MS2_Data)):
                if self.MS2_Data[i]['Data'].MS2_Data['Pre_MZ'][ii] == self.MS2_Data[i]['Data'].MS2_Data['Pre_MZ'][ii-1]:
                    V_Place = list(filter(lambda x:V_List[x]==self.MS2_Data[i]['Data'].MS2_Data['CE'][ii-1],range(len(V_List))))[0]
                    if V_Place == len(V_List)-1:
                        self.MS2_Data[i]['Data'].MS2_Data['CE'][ii] = V_List[0]
                    else:
                        self.MS2_Data[i]['Data'].MS2_Data['CE'][ii] = V_List[V_Place+1]
                else:
                    self.MS2_Data[i]['Data'].MS2_Data['CE'][ii] = V_List[0]
            merge_MS2_Data = pd.DataFrame(columns=['Scan_Time','Pre_MZ','MZ_List','Int_List','Rel_Int','CE','DataName','Pre-Int'])
            temp_MS2_Data = pd.DataFrame(columns=['Scan_Time','Pre_MZ','MZ_List','Int_List','Rel_Int','CE','DataName','Pre-Int'])
            for ii in range(1,len(self.MS2_Data[i]['Data'].MS2_Data)):
                if len(temp_MS2_Data) == 0:
                    temp_MS2_Data = pd.concat([temp_MS2_Data,self.MS2_Data[i]['Data'].MS2_Data.iloc[[ii],:]])
                    temp_MS2_Data.reset_index(drop=True,inplace=True)
                else:
                    CE_place = list(filter(lambda x:V_List[x] == self.MS2_Data[i]['Data'].MS2_Data['CE'][ii],range(len(V_List))))
                    if CE_place[0]<len(V_List)-1:
                        temp_MS2_Data = pd.concat([temp_MS2_Data,self.MS2_Data[i]['Data'].MS2_Data.iloc[[ii],:]])
                        temp_MS2_Data.reset_index(drop=True,inplace=True)
                    else:
                        temp_MS2_Data = pd.concat([temp_MS2_Data,self.MS2_Data[i]['Data'].MS2_Data.iloc[[ii],:]])
                        temp_MS2_Data.reset_index(drop=True,inplace=True)
                        MS2_Int_temp = []
                        MS2_MZ_temp = []
                        Scan_Time = []
                        '''   '''
                        for iii in range(len(temp_MS2_Data)):
                            Scan_Time.append(temp_MS2_Data['Scan_Time'][iii])
                            if len(MS2_Int_temp) == 0:
                                for iv in range(len(temp_MS2_Data['MZ_List'][iii])):
                                    MS2_Int_temp.append(temp_MS2_Data['Int_List'][iii][iv])
                                    MS2_MZ_temp.append(temp_MS2_Data['MZ_List'][iii][iv])
                            else:
                                for iv in range(len(temp_MS2_Data['MZ_List'][iii])):
                                    MS2_match = list(filter(lambda x:abs(temp_MS2_Data['MZ_List'][iii][iv]-MS2_MZ_temp[x])/temp_MS2_Data['MZ_List'][iii][iv]<float(self.LineEdit_MS_Tor.text())/1000000,range(len(MS2_MZ_temp))))
                                    if len(MS2_match) > 0:
                                        MS2_Int_temp[MS2_match[0]] = max(MS2_Int_temp[MS2_match[0]],temp_MS2_Data['Int_List'][iii][iv])
                                    else:
                                        MS2_Int_temp.append(temp_MS2_Data['Int_List'][iii][iv])
                                        MS2_MZ_temp.append(temp_MS2_Data['MZ_List'][iii][iv])
                        temp_range = pd.DataFrame({'MZ_List':MS2_MZ_temp,'Int_List':np.array(MS2_Int_temp)/max(MS2_Int_temp)})
                        temp_range.sort_values(by='MZ_List',inplace=True)
                        temp_range.reset_index(drop=True,inplace=True)
                        MS2_MZ_temp = list(temp_range['MZ_List'])
                        MS2_Int_temp = list(temp_range['Int_List'])
                        temp_merge_MS2_Data = pd.DataFrame({'Scan_Time':[np.mean(Scan_Time)],'Pre_MZ':[temp_MS2_Data['Pre_MZ'][0]],'MZ_List':[MS2_MZ_temp],'Int_List':[np.array(MS2_Int_temp)/max(MS2_Int_temp)],'Rel_Int':[max(MS2_Int_temp)],'CE':['merge'],'DataName':[temp_MS2_Data['DataName'][0]],'Pre-Int':[temp_MS2_Data['Pre-Int'][0]]})
                        merge_MS2_Data = pd.concat([merge_MS2_Data,temp_merge_MS2_Data])
                        merge_MS2_Data.reset_index(drop=True,inplace=True)
                        temp_MS2_Data = pd.DataFrame(columns=['Scan_Time','Pre_MZ','MZ_List','Int_List','Rel_Int','CE','DataName','Pre-Int'])
            self.MS2_Data[i]['Data'].MS2_Data = merge_MS2_Data.copy()
        self.Align.RefList['MS2_MZ'] = self.Align.RefList['MS2_MZ'].astype('object')
        self.Align.RefList['MS2_Int'] = self.Align.RefList['MS2_Int'].astype('object')
        self.Label_process_sub.setText('MS2 assine')
        self.process_bar.setMaximum(len(self.Align.RefList))
        self.process_bar.setValue(0)
        QApplication.processEvents()
        self.MergeMS2 = []
        for i in self.MS2_Data.keys():
            self.MS2_Data[i]['Data'].Calculate_MS2_RI()
            if len(self.MergeMS2) == 0:
                self.MergeMS2 = self.MS2_Data[i]['Data'].MS2_Data.copy()
            else:
                self.MergeMS2 = pd.concat([self.MergeMS2,self.MS2_Data[i]['Data'].MS2_Data])
        if len(self.MergeMS2)>0:
            self.MergeMS2.sort_values(['Pre_MZ','Scan_RI'],ascending=(True),inplace=True)
            self.MergeMS2.reset_index(drop=True,inplace=True)
        temp_MZ_Tor = self.Align.AlignmentParam['MZ_Tor']
        for i in range(len(self.Align.RefList)):
            self.process_bar.setValue(self.process_bar.value()+1)
            QApplication.processEvents()
            Valified_List = list(filter(lambda x:abs(self.Align.RefList.at[i,'m/z']-self.MergeMS2.at[x,'Pre_MZ'])/self.Align.RefList.at[i,'m/z']<temp_MZ_Tor and self.Align.RefList.at[i,'RI']*(1-float(self.LineEdit_RT_Tor.text())/100)<self.MergeMS2.at[x,'Scan_RI']<self.Align.RefList.at[i,'RI']*((1+float(self.LineEdit_RT_Tor.text())/100)),range(len(self.MergeMS2))))
            if len(Valified_List) >0:
                Valified_Index = np.where(abs(np.array(self.MergeMS2.loc[Valified_List,'Scan_RI'])-self.Align.RefList.at[i,'RI'])==min(abs(np.array(self.MergeMS2.loc[Valified_List,'Scan_RI'])-self.Align.RefList.at[i,'RI'])))[0][0]
                self.Align.RefList.at[i,'MS2_MZ'] = self.MergeMS2.at[Valified_List[Valified_Index],'MZ_List'].copy()
                self.Align.RefList.at[i,'MS2_Int'] = self.MergeMS2.at[Valified_List[Valified_Index],'Int_List'].copy()

        ''' Filter '''
        for i in ClassList:
            self.Align.RefList['Group '+str(i)] = 0
            for ii in range(len(self.Align.RefList)):
                Number_Class = 0
                for iii in self.SampleClass.keys():
                    if self.SampleClass[iii] == i:
                        Number_Class += 1
                        if self.Align.RefList[iii][ii] >0:
                            self.Align.RefList['Group '+str(i)][ii] += 1
                self.Align.RefList['Group '+str(i)][ii] = self.Align.RefList['Group '+str(i)][ii]/Number_Class
        #self.Align.RefList.to_excel(self.filepath_title+'Alignment-Origin-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        self.Align.BlankFilter()
        self.Align.RefList.to_excel(self.filepath_title+'Alignment-BlankRemove-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        #self.Align.RefList.reset_index(drop=True,inplace=True)
        #self.Align.RefList.to_excel(self.Align_path,index=False)
        mgf_output = ''
        csv_mgf = pd.DataFrame(columns=['row ID','row m/z','row retention time','correlation group ID','annotation network number','best ion','auto MS2 verify','identified by n=','partners','neutral M mass','Peak height'])
        for i in  range(len(self.Align.RefList)):
            MZ = self.Align.RefList['m/z'][i]
            RT = np.around(self.Align.RefList['RT'][i]/60,3)
            Int = self.Align.RefList['Int'][i]
            temp_csv_mgf = pd.DataFrame([(i+1,MZ,RT,'','','','','','','',Int)],columns=['row ID','row m/z','row retention time','correlation group ID','annotation network number','best ion','auto MS2 verify','identified by n=','partners','neutral M mass','Peak height'])
            csv_mgf = pd.concat([csv_mgf,temp_csv_mgf])
            if self.Polarity == 'Positive':
                mgf_output = mgf_output+'BEGIN IONS\nFEATURE_ID='+str(i+1)+'\nPEPMASS='+str(MZ)+'\nSCANS='+str(i+1)+'\nRTINSECONDS='+str(RT)+'\nCHARGE=1+\nMSLEVEL=2\n'
            elif self.Polarity == 'Negative':
                mgf_output = mgf_output+'BEGIN IONS\nFEATURE_ID='+str(i+1)+'\nPEPMASS='+str(MZ)+'\nSCANS='+str(i+1)+'\nRTINSECONDS='+str(RT)+'\nCHARGE=1-\nMSLEVEL=2\n'
            if type(self.Align.RefList['MS2_MZ'][i]) == np.ndarray:
                for ii in range(len(self.Align.RefList['MS2_MZ'][i])):
                    mgf_output = mgf_output+str(self.Align.RefList['MS2_MZ'][i][ii])+' '+str(self.Align.RefList['MS2_Int'][i][ii])+'\n'
            mgf_output = mgf_output+'END IONS\n\n'
        with open(self.filepath_title+'Alignment-MS2-MGF-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.mgf','w')as mgfFile:
            mgfFile.write(mgf_output)
        csv_mgf.to_csv(self.filepath_title+'Alignment-MS2-CSV-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.csv',index=False)
        self.Label_process_sub.setText('Finished')
        self.process_bar.setMaximum(100)
        self.process_bar.setValue(100)
        
        
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
        
class EazyMZDataProcess(object):   
    def __init__(self,DataPath):
        if DataPath.rfind('/') != -1 :
            self.DataName = DataPath[DataPath.rfind('/')+1:len(DataPath)-5]
        else:
            self.DataName = DataPath[DataPath.rfind('\\')+1:len(DataPath)-5]
        self.OriginData = pyopenms.MSExperiment()
        self.file_path = DataPath
        ''' 储存文件pyopenms.MzMLFile().store("filtered.mzML", exp) '''
        if self.file_path.endswith('mzML'):
            pyopenms.MzMLFile().load(self.file_path,self.OriginData)
        elif self.file_path.endswith('mzXML'):
            pyopenms.MzXMLFile().load(self.file_path,self.OriginData)
        self.OriginData.sortSpectra(True)
        self.__param = {'MS1_Tor':0.000010,'RT_Tor':2,'min_Int':10000,'min_RT':6,'max_Noise':2000,
                        'Deconvolution':False,'FeatureDetectPlot':3,'MergeRule':'Intersection',
                        'UpDown_gap':10,'saveAutoList':False,'smooth':5,'Points':25}
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
                temp_range = list(filter(lambda x:MZ_temp[x]<=Pre_temp,range(len(MZ_temp)))) # and Int_temp[x]>1000
                temp_range = list(filter(lambda x:len(np.where(abs(MZ_temp[temp_range]-MZ_temp[x])/MZ_temp[x]<self.__param['MS1_Tor']))==1 or
                                         Int_temp[x]==max(Int_temp[np.where(abs(MZ_temp[temp_range]-MZ_temp[x])/MZ_temp[x]<self.__param['MS1_Tor'])]),temp_range))
                if len(temp_range)>0:
                    MZ_temp = MZ_temp[temp_range]
                    Int_temp = Int_temp[temp_range]
                    self.MS2_Pre.append(Pre_temp)
                    self.MS2_RT_List.append(i.getRT())
                    self.MS2_MZ_List.append(np.around(MZ_temp,5))
                    self.MS2_Int_List.append(np.around(Int_temp,0))
                    self.MS2_RelInt.append(max(Int_temp[temp_range]))
                else:
                    MZ_temp = []
                    Int_temp = []
                    self.MS2_Pre.append(Pre_temp)
                    self.MS2_RT_List.append(i.getRT())
                    self.MS2_MZ_List.append(np.around(MZ_temp,5))
                    self.MS2_Int_List.append(np.around(Int_temp,0))
                    self.MS2_RelInt.append(0)
        self.Origin_RT_List = np.around(self.Origin_RT_List,3)
        self.MS2_RT_List = np.around(self.MS2_RT_List,3)
        self.MS1_Data = {'Scan_Time':self.Origin_RT_List,'MZ_List':self.Origin_MZ_List,'Int_List':self.Origin_Int_List}
        self.MS1_Data = pd.DataFrame(self.MS1_Data)
        self.MS2_Data = {'Scan_Time':self.MS2_RT_List,'Pre_MZ':self.MS2_Pre,'MZ_List':self.MS2_MZ_List,'Int_List':self.MS2_Int_List,'Rel_Int':self.MS2_RelInt}
        self.MS2_Data = pd.DataFrame(self.MS2_Data)
        self.__param['Flow_RT'] = int(self.__param['min_RT']/(4*sum(np.diff(self.Origin_RT_List))/len(np.diff(self.Origin_RT_List))))+1
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
        Heuristic_peak_list['RI_left'] = Heuristic_peak_list['RI'].apply(lambda x:x*0.98)
        Heuristic_peak_list['RI_right'] = Heuristic_peak_list['RI'].apply(lambda x:x*1.02)
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
                        if abs(i_Ref_MZ-i_Sample_MZ)/i_Ref_MZ<self.__param['MS1_Tor'] and abs(i_Ref_Time-i_Sample_Time)/i_Ref_Time<self.__param['RT_Tor']:
                            Score_matrix[i_Ref,i_Sample] = 0.5*np.exp(-0.5*((i_Sample_Time-i_Ref_Time)/(i_Ref_Time*self.__param['RT_Tor']))**2)+(1-0.5)*np.exp(-0.5*((i_Sample_MZ-i_Ref_MZ)/(i_Sample_MZ*self.__param['MS1_Tor']))**2)
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
        Alltemp_HPL['RI_left'] = Alltemp_HPL['RI'].apply(lambda x:x*(1-self.__param['RT_Tor']))
        Alltemp_HPL['RI_right'] = Alltemp_HPL['RI'].apply(lambda x:x*(1+self.__param['RT_Tor']))
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
                    MZ_place = EazyMZDataProcess.ClosestPosition(MZ,self.Origin_MZ_List[i])
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
                    if not temp_int:
                        Int_List.append(0)
                        #MZ_List.append(0)
                    else:
                        temp_int = np.array(temp_int)
                        Int_List.append(temp_int.max())
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
        if plt_for_test == True:
            plt.figure()
            plt.title(str(MZ))
            plt.plot(RT_List,Int_List)
            plt.figure()
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
                    if not temp_int:
                        Int_List.append(0)
                    else:
                        temp_int = np.array(temp_int)
                        Int_List.append(temp_int.max())
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
        if plt_for_test == True:
            plt.figure()
            plt.title(str(MZ))
            plt.plot(RT_List,Int_List)
            plt.figure()
        return RT_List, Int_List
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
        n_place = np.where(RIIS['C']<=RI/100)[0]
        n1_place = np.where(RIIS['C']>RI/100)[0]
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
        RT = (RI/100-C_number)*(t_n1-t_n)+t_n
        return RT
    def RT_to_RI(RT,RIIS):
        n_place = np.where(RIIS['RT']<=RT)[0]
        n1_place = np.where(RIIS['RT']>RT)[0]
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
        RI = 100*(C_number+(RT-t_n)/(t_n1-t_n))
        return RI
    def Heuristic_PeakDetect(self):
        warnings.filterwarnings('ignore')
        def RI_to_RT(RI,RIIS):
            n_place = np.where(RIIS['C']<=RI/100)[0]
            n1_place = np.where(RIIS['C']>RI/100)[0]
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
            RT = (RI/100-C_number)*(t_n1-t_n)+t_n
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
            n_place = np.where(RIIS['RT']<=RT)[0]
            n1_place = np.where(RIIS['RT']>RT)[0]
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
            RI = 100*(C_number+(RT-t_n)/(t_n1-t_n))
            return RI
        Heuristic_peak_list = self.Heuristic_peak_list.copy()
        Final_Peak_Detect = pd.DataFrame(columns=['AverageMZ', 'RT', 'RI', 'Int','RTList','SimilarityScore'])
        Heuristic_peak_list['ActualNumber'] = 0
        for ii in range(len(Heuristic_peak_list)):   
            MZ = Heuristic_peak_list.at[ii,'AverageMZ']
            RI_left = Heuristic_peak_list.at[ii,'RI_left']
            RI_right = Heuristic_peak_list.at[ii,'RI_right']
            RI_List = Heuristic_peak_list.at[ii,'RI_List']
            RTL = RI_to_RT(RI_left,self.RIIS)
            RTR = RI_to_RT(RI_right,self.RIIS)
            #Gradient_Code = Heuristic_peak_list.at[ii,'Gradient Code']
            [RT_List,Int_List] = self.ExtractDataPoint(MZ,self.Origin_RT_List[-1],0,smooth_index=0)
            Noise = sum(Int_List)/len(Int_List)+1
            [RT_List,Int_List_Origin] = self.ExtractDataPoint(MZ,RTR+10,RTL-10,smooth_index=0)
            [RT_List,Int_List] = self.ExtractDataPoint(MZ,RTR+10,RTL-10)
            Int_List = list(map(lambda x:sum(Int_List[x-2:x+2])/len(Int_List[x-2:x+2]) if x in range(2,len(Int_List)-2) else Int_List[x],range(len(Int_List))))
            if max(Int_List)>0 and max(Int_List_Origin)>0:
                Int_Correction = max(Int_List_Origin)/max(Int_List)
            else:
                Int_Correction = 1
            FD_List = np.array(list(map(lambda x: EazyMZDataProcess.TFFD(x,Int_List,RT_List), range(2, len(Int_List)-2))))
            FD_List = np.array(list(map(lambda x:GaussSmooth(FD_List[x-2:x+2+1]) if x in range(2,len(FD_List)-2) else FD_List[x],range(len(FD_List)))))
            SD_List = np.array(list(map(lambda x: EazyMZDataProcess.TFSD(x,Int_List,RT_List), range(2, len(Int_List)-2))))
            SD_List = np.array(list(map(lambda x:GaussSmooth(SD_List[x-2:x+2+1]) if x in range(2,len(SD_List)-2) else SD_List[x],range(len(SD_List)))))
            ABS_FD_List = abs(np.array(FD_List))
            if Heuristic_peak_list.at[ii,'PeakNumber'] == 1:
                FD_Median = EazyMZDataProcess.FD_Line(ABS_FD_List)
                FD_Median_N = FD_Median*(-1)
                SD_Median = min(SD_List)*0.1
            else:
                FD_Median = EazyMZDataProcess.FD_Line(ABS_FD_List)
                FD_Median = FD_Median*Heuristic_peak_list.at[ii,'PeakNumber']*Heuristic_peak_list.at[ii,'PeakNumber']*self.__param['min_Int']/max(Int_List)
                FD_Median_N = FD_Median*(-1)
                SD_Median = sum(SD_List[(min(SD_List)*0.1<=SD_List) & (SD_List<=0)])/max(len(SD_List[(min(SD_List)*0.1<=SD_List) & (SD_List<=0)]),1)
                #SD_Median = min(SD_List)*0.1*Heuristic_peak_list.at[ii,'PeakNumber']*self.__param['min_Int']/max(Int_List)
            FD_P = list(filter(lambda x: FD_List[x] > FD_Median, range(1, len(FD_List))))
            FD_N = list(filter(lambda x: FD_List[x] < FD_Median_N, range(1, len(FD_List))))
            FD_Change = list(filter(lambda x: FD_List[x-1] > 0 and FD_List[x] < 0, range(1, len(FD_List))))
            SD_Place = list(filter(lambda x: SD_List[x] < SD_Median, range(len(SD_List))))
            Begin_Place,Complet_BP = EazyMZDataProcess.Find_FContinuous(FD_P, FD_P, FDMode='P',mergeRule='Intersection',FC_Number=2)
            Begin_Place += 2
            Complet_BP += 2
            End_Place,Complet_EP = EazyMZDataProcess.Find_FContinuous(FD_N, FD_N, FDMode='N',mergeRule='Intersection',FC_Number=2)
            End_Place += 2
            Complet_EP += 2
            Peak_Place = EazyMZDataProcess.Find_SDChange(FD_Change, SD_Place)
            Peak_Place += 2
            FD_Sequence = np.concatenate((Begin_Place, End_Place))
            FD_Sequence.sort()
            if len(Begin_Place) >= 1 and len(End_Place) >= 1 and Begin_Place[0] < End_Place[-1]:
                PeakDetect_Table = pd.DataFrame({'Begin_Place':Begin_Place})
                PeakDetect_Table['Peak_Place'] = PeakDetect_Table['Begin_Place'].apply(lambda x:[])
                PeakDetect_Table['End_Place'] = PeakDetect_Table['Begin_Place'].apply(lambda x:[])
                PeakDetect_Table['End_Place'] = PeakDetect_Table['End_Place'].astype('object')
                PeakDetect_Table['Peak_Place'] = PeakDetect_Table['Peak_Place'].astype('object')
                for i_FD_Seq in range(len(Begin_Place)):
                    if i_FD_Seq < len(Begin_Place)-1:
                        Valified_Peak_Place = list(filter(lambda x:PeakDetect_Table.at[i_FD_Seq,'Begin_Place']<Peak_Place[x]<PeakDetect_Table.at[i_FD_Seq+1,'Begin_Place'],range(len(Peak_Place))))
                        if len(Valified_Peak_Place) == 0:
                            continue
                        PeakDetect_Table.at[i_FD_Seq,'Peak_Place'] = list(Peak_Place[Valified_Peak_Place])
                        Valified_End_Place = list(filter(lambda x:PeakDetect_Table.at[i_FD_Seq,'Peak_Place'][0]<End_Place[x]<=PeakDetect_Table.at[i_FD_Seq+1,'Begin_Place'],range(len(End_Place))))
                        PeakDetect_Table.at[i_FD_Seq,'End_Place'] = list(End_Place[Valified_End_Place])
                    else:
                        Valified_Peak_Place = list(filter(lambda x:PeakDetect_Table.at[i_FD_Seq,'Begin_Place']<Peak_Place[x],range(len(Peak_Place))))
                        if len(Valified_Peak_Place) == 0:
                            continue
                        PeakDetect_Table.at[i_FD_Seq,'Peak_Place'] = list(Peak_Place[Valified_Peak_Place])
                        Valified_End_Place = list(filter(lambda x:End_Place[x]>PeakDetect_Table.at[i_FD_Seq,'Peak_Place'][0],range(len(End_Place))))
                        PeakDetect_Table.at[i_FD_Seq,'End_Place'] = list(End_Place[Valified_End_Place])
                Valified = []
                for i_FD_Seq in range(len(PeakDetect_Table)):
                    if len(PeakDetect_Table.at[i_FD_Seq,'Peak_Place'])==0 or len(PeakDetect_Table.at[i_FD_Seq,'End_Place'])==0:
                        continue
                    Valified.append(i_FD_Seq)
                    ii_FD_VPP = []
                    ii_FD_VEP = []
                    for ii_FD_Seq in range(len(PeakDetect_Table.at[i_FD_Seq,'End_Place'])):
                        Valified_Peak_Place = list(filter(lambda x:PeakDetect_Table.at[i_FD_Seq,'End_Place'][ii_FD_Seq]>PeakDetect_Table.at[i_FD_Seq,'Peak_Place'][x],range(len(PeakDetect_Table.at[i_FD_Seq,'Peak_Place']))))
                        Valified_Peak_Place = list(filter(lambda x:Int_List[Peak_Place[Valified_Peak_Place[x]]]==max(np.array(Int_List)[Peak_Place[Valified_Peak_Place]]),range(len(Valified_Peak_Place))))[0]
                        ii_FD_VPP.append(Valified_Peak_Place)
                        ii_FD_VEP.append(abs(1-(PeakDetect_Table.at[i_FD_Seq,'Peak_Place'][Valified_Peak_Place]-PeakDetect_Table.at[i_FD_Seq,'Begin_Place'])/(PeakDetect_Table.at[i_FD_Seq,'End_Place'][ii_FD_Seq]-PeakDetect_Table.at[i_FD_Seq,'Peak_Place'][Valified_Peak_Place])))
                    Final_Valified_Place = list(filter(lambda x:ii_FD_VEP[x]==min(ii_FD_VEP),range(len(ii_FD_VEP))))[0]
                    PeakDetect_Table.at[i_FD_Seq,'Peak_Place'] = PeakDetect_Table.at[i_FD_Seq,'Peak_Place'][ii_FD_VPP[Final_Valified_Place]]
                    PeakDetect_Table.at[i_FD_Seq,'End_Place'] = PeakDetect_Table.at[i_FD_Seq,'End_Place'][Final_Valified_Place]
                PeakDetect_Table = PeakDetect_Table.iloc[Valified,:].copy()
                PeakDetect_Table.reset_index(drop=True,inplace=True)
                for i_wave in range(len(PeakDetect_Table)):
                    if PeakDetect_Table.at[i_wave,'Begin_Place'] >= int(self.__param['Points']/3):
                        Begin_min = list(filter(lambda x:min(Int_List[PeakDetect_Table.at[i_wave,'Begin_Place']-int(self.__param['Points']/3):PeakDetect_Table.at[i_wave,'Begin_Place']])==Int_List[x],range(PeakDetect_Table.at[i_wave,'Begin_Place']-int(self.__param['Points']/3),PeakDetect_Table.at[i_wave,'Begin_Place'])))
                        PeakDetect_Table.at[i_wave,'Begin_Place'] = Begin_min[0]
                    if PeakDetect_Table.at[i_wave,'End_Place'] < len(Int_List)-int(self.__param['Points']/3):
                        End_min = list(filter(lambda x:min(Int_List[PeakDetect_Table.at[i_wave,'End_Place']:PeakDetect_Table.at[i_wave,'End_Place']+int(self.__param['Points']/3)])==Int_List[x],range(PeakDetect_Table.at[i_wave,'End_Place'],PeakDetect_Table.at[i_wave,'End_Place']+int(self.__param['Points']/3))))
                        PeakDetect_Table.at[i_wave,'End_Place'] = End_min[0]
                Peak_Begin = list(PeakDetect_Table['Begin_Place'])
                Peak_End = list(PeakDetect_Table['End_Place'])
                Peak_Top = list(PeakDetect_Table['Peak_Place'])
                temp_Final_Peak = pd.DataFrame(columns=['AverageMZ', 'RT', 'RI', 'Int','RTList','SimilarityScore'])
                for v in range(len(Peak_Begin)):
                    RT = RT_List[Peak_Top[v]]
                    RI = RT_to_RI(RT, self.RIIS)
                    Int = max(Int_List_Origin[Peak_Top[v]],Int_List[Peak_Top[v]]*Int_Correction)
                    #[BaseValue,Noise] = self.Baseline(MZ,RT_List[Peak_End[v]+1],RT_List[Peak_Begin[v]])
                    RTList = RT_List[Peak_Begin[v]:Peak_End[v]+1]
                    Peak_Int_List = Int_List[Peak_Begin[v]:Peak_End[v]+1]
                    #rate = max(Peak_Int_List[0],Peak_Int_List[-1])/Int
                    if Int > self.__param['min_Int'] and Int/Noise>=3 and len(RTList)>=self.__param['Points'] and len(list(filter(lambda x:(1-self.__param['RT_Tor'])<RI/x<(1+self.__param['RT_Tor']),RI_List)))>=1: 
                        #IntList = Int_List[Peak_Begin[v]:Peak_End[v]+1]
                        RTPlace = list(filter(lambda x:Peak_Int_List[x]==max(Peak_Int_List),range(len(Peak_Int_List))))
                        Int_List_Gaussian_1 = list(map(lambda x:gaussian(x,b=RTList[RTPlace[0]],c=1),RTList))
                        Int_List_Gaussian_2 = list(map(lambda x:gaussian(x,b=RTList[RTPlace[0]],c=1.3),RTList))
                        Int_List_Gaussian_3 = list(map(lambda x:gaussian(x,b=RTList[RTPlace[0]],c=1.6),RTList))
                        Int_List_Gaussian_4 = list(map(lambda x:gaussian(x,b=RTList[RTPlace[0]],c=2),RTList))
                        SimilarityScore = max(EazyMZDataProcess.CosineSimilarity(Peak_Int_List,Int_List_Gaussian_1),
                                         EazyMZDataProcess.CosineSimilarity(Peak_Int_List,Int_List_Gaussian_2),
                                         EazyMZDataProcess.CosineSimilarity(Peak_Int_List,Int_List_Gaussian_3),
                                         EazyMZDataProcess.CosineSimilarity(Peak_Int_List,Int_List_Gaussian_4))
                        Same_list = list(filter(lambda x:abs(Final_Peak_Detect['AverageMZ'][x]-MZ)/MZ<self.__param['MS1_Tor'] and abs(Final_Peak_Detect['RT'][x]-RT)<2,range(len(Final_Peak_Detect))))
                        if len(Same_list) == 0:
                            temp_Final_Peak_v = pd.DataFrame([(MZ, RT, RI, Int,RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'RI', 'Int','RTList','SimilarityScore'])
                            temp_Final_Peak = pd.concat([temp_Final_Peak,temp_Final_Peak_v])
                            temp_Final_Peak.reset_index(drop=True,inplace=True)
                Score_matrix= np.zeros([len(RI_List),len(temp_Final_Peak)])
                for i_Sm_row in range(len(RI_List)):
                    for i_Sm_col in range(len(temp_Final_Peak)):
                        if abs(temp_Final_Peak.at[i_Sm_col,'RI']/RI_List[i_Sm_row]-1) <=0.02:
                            Score_matrix[i_Sm_row,i_Sm_col] = temp_Final_Peak.at[i_Sm_col,'Int']
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
                    temp_Final_Peak.sort_values('Int',inplace=True)
                    temp_Final_Peak.reset_index(drop=True,inplace=True)
                    temp_Final_Peak = temp_Final_Peak.iloc[range(Heuristic_peak_list.at[ii,'PeakNumber']),:].copy()
                Final_Peak_Detect = pd.concat([Final_Peak_Detect, temp_Final_Peak])
                Final_Peak_Detect.reset_index(drop=True,inplace=True)
        return [Final_Peak_Detect.copy(),Heuristic_peak_list.copy()]

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
        ABS_FD_List = list(filter(lambda x: x < max(ABS_FD_List)*0.10 and x > 0, ABS_FD_List))
        FD_Median = np.median(ABS_FD_List)
        return FD_Median
    def SD_Line(SD_List):
        SD_List = abs(SD_List)
        if len(SD_List) >= 1:
            SD_Median = min(SD_List)*0.05
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
    def Find_FContinuous(FD_P, Diff_P, FDMode='P' ,FC_Number=2 ,mergeRule='Intersection'):
        if mergeRule == 'Intersection':
            First_List = list(filter(lambda x:x in FD_P,Diff_P))
        elif mergeRule == 'Union':
            First_List = list(set(FD_P+Diff_P))
            First_List.sort()
        if len(First_List)>0:
            Second_List = [First_List[0]]
            Second_Score = [First_List[0]]
            for i_FSB in range(1,len(First_List)):
                if First_List[i_FSB]-First_List[i_FSB-1]==1:
                    Second_Score[-1] = First_List[i_FSB]
                else:
                    Second_List.append(First_List[i_FSB])
                    Second_Score.append(First_List[i_FSB])
            Second_List = np.array(Second_List)
            Second_Score = np.array(Second_Score)
            temp_Place = np.where(Second_Score-Second_List >= FC_Number)[0]
            if len(temp_Place) > 0:
                if FDMode == 'P':
                    return Second_List[temp_Place],Second_Score[temp_Place]
                elif FDMode == 'N':
                    return Second_Score[temp_Place],Second_List[temp_Place]
            else:
                if len(FD_P) <= FC_Number+1 or len(Diff_P) <= FC_Number+1:
                    Second_List_F = [FD_P[0]]
                    Second_Score_F = [FD_P[0]]
                    for i_FSB in range(1,len(FD_P)):
                        if FD_P[i_FSB]-FD_P[i_FSB-1]==1:
                            Second_Score_F[-1] = FD_P[i_FSB]
                        else:
                            Second_List_F.append(FD_P[i_FSB])
                            Second_Score_F.append(FD_P[i_FSB])
                    Second_List_F = np.array(Second_List_F)
                    Second_Score_F = np.array(Second_Score_F)
                    temp_Place_F = np.where(Second_Score_F-Second_List_F >= FC_Number)[0]
                    Second_List_D = [Diff_P[0]]
                    Second_Score_D = [Diff_P[0]]
                    for i_FSB in range(1,len(Diff_P)):
                        if Diff_P[i_FSB]-Diff_P[i_FSB-1]==1:
                            Second_Score_D[-1] = Diff_P[i_FSB]
                        else:
                            Second_List_D.append(Diff_P[i_FSB])
                            Second_Score_D.append(Diff_P[i_FSB])
                    Second_List_D = np.array(Second_List_D)
                    Second_Score_D = np.array(Second_Score_D)
                    temp_Place_D = np.where(Second_Score_D-Second_List_D >= FC_Number)[0]
                    if len(temp_Place_F) > 0 and len(temp_Place_D) > 0:
                        if FDMode == 'P':
                            return min(Second_List_F[temp_Place_F],Second_List_D[temp_Place_D]),max(Second_Score_F[temp_Place_F],Second_Score_D[temp_Place_D])
                        elif FDMode == 'N':
                            return max(Second_Score_F[temp_Place_F],Second_Score_D[temp_Place_D]),min(Second_List_F[temp_Place_F],Second_List_D[temp_Place_D])
                    else:  
                        return np.array([]),np.array([])   
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
        self.Final_Peak_Detect['S/N']=self.Final_Peak_Detect['AverageMZ'].apply(lambda x:0)
        for i in range(len(self.Final_Peak_Detect)):
            [RT_List,Int_List] = EazyMZDataProcess.ExtractDataPoint(self,self.Final_Peak_Detect['AverageMZ'][i],self.Origin_RT_List[-1],0,smooth_index=0)
            Noise = sum(np.array(Int_List)[Int_List<max(Int_List)*0.8])/len(np.array(Int_List)[Int_List<max(Int_List)*0.8])+1
            self.Final_Peak_Detect.loc[i,'S/N'] = self.Final_Peak_Detect['Int'][i]/Noise
        if drop == True:
            self.Final_Peak_Detect = self.Final_Peak_Detect[self.Final_Peak_Detect['S/N']>Threshold].copy()
            self.Final_Peak_Detect.reset_index(drop=True,inplace=True)
    def set_RIIS(self,RIIS):
        self.RIIS = pd.read_excel(RIIS)
        if 'RT' not in list(self.RIIS.keys()):
            self.RIIS['RT']=0.00
            RTL = 60
            RTR = self.Origin_RT_List[-1]
            for ii in range(len(self.RIIS)):                
                MZ = self.RIIS['m/z'][ii]
                [RT_List, Int_List]=self.ExtractDataPoint(MZ,RTR,RTL)
                RT_place = np.where(Int_List==max(Int_List))[0]
                self.RIIS['RT'][ii] = RT_List[RT_place[0]].copy()
            self.RIIS['RT']=0.00
            self.RIIS['Candidate_RT'] = ''
            self.RIIS['Candidate_Int'] = ''
            RTL = 0
            RTR = self.Origin_RT_List[-1]
            for i in range(len(self.RIIS)):                
                MZ = self.RIIS['m/z'][i]
                [RT_List, Int_List]=self.ExtractDataPoint(MZ,RTR,RTL)
                RT_place_1 = np.where(Int_List==max(Int_List))[0][0]
                Int_List_d = Int_List[0:RT_place_1-int(self.__param['Points']/2)]+Int_List[RT_place_1+int(self.__param['Points']/2):-1]
                RT_place_2 = np.where(Int_List_d==max(Int_List_d))[0][0]
                if RT_place_2 >= RT_place_1-int(self.__param['Points']/2):
                    RT_place_2 += int(self.__param['Points']/2)*2
                self.RIIS['Candidate_RT'][i] = [RT_List[RT_place_1],RT_List[RT_place_2]]
                self.RIIS['Candidate_Int'][i] = [Int_List[RT_place_1],Int_List[RT_place_2]]
                self.RIIS['RT'][i] = RT_List[RT_place_1].copy()
            for i in range(len(self.RIIS.iloc[:,0])):
                if i == 0:
                    if self.RIIS['RT'][i] > self.RIIS['RT'][i+1]:
                        self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][1]
                        if self.RIIS['RT'][i]>self.RIIS['RT'][i+1]:
                            self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][0]
                elif i == len(self.RIIS.iloc[:,0])-1:
                    if self.RIIS['RT'][i] < self.RIIS['RT'][i-1]:
                        self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][1]
                        if self.RIIS['RT'][i]<self.RIIS['RT'][i-1]:
                            self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][0]
                else:
                    if self.RIIS['RT'][i] < self.RIIS['RT'][i-1] or self.RIIS['RT'][i]>self.RIIS['RT'][i+1]:
                        self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][1]
                        if self.RIIS['RT'][i]<self.RIIS['RT'][i-1] or self.RIIS['RT'][i]>self.RIIS['RT'][i+1]:
                            self.RIIS['RT'][i] = self.RIIS['Candidate_RT'][i][0]
            for i in range(1,len(self.RIIS.iloc[:,1])-1):
                if min(self.RIIS['Candidate_Int'][i])/max(self.RIIS['Candidate_Int'][i])>0.4 and self.RIIS['RT'][i] != self.RIIS['Candidate_RT'][i][1]:
                    if self.RIIS.at[i-1,'RT']<self.RIIS['Candidate_RT'][i][1]<self.RIIS.at[i+1,'RT']:
                        score_1 = (self.RIIS['Candidate_RT'][i][0]-self.RIIS.at[i-1,'RT'])/(self.RIIS.at[i+1,'RT']-self.RIIS['Candidate_RT'][i][0])
                        score_2 = (self.RIIS['Candidate_RT'][i][1]-self.RIIS.at[i-1,'RT'])/(self.RIIS.at[i+1,'RT']-self.RIIS['Candidate_RT'][i][1])
                        if abs(1-score_2)<abs(1-score_1):
                            self.RIIS.at[i,'RT']=self.RIIS.at[i,'Candidate_RT'][1]
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

class DataAlignment(object):
    def __init__(self,Label_process_sub,process_bar):
        self.Label_process_sub = Label_process_sub
        self.process_bar = process_bar
        self.DataBase = pd.DataFrame(columns=['Data','Data_Name','Tag'])
        self.AlignmentParam = {'MZ_Tor':0.000010,'RT_Tor':0.02,'A':0.75,
                               'Miss_Filter':0.8,'Threshold':5,'RI_Alignment':False,
                               'RT_min_Tor':6}
        self.RefList=pd.DataFrame(columns=['m/z','RT','MS_List','RT_List'])
        self.RefList['m/z'] = self.RefList['m/z'].map(lambda x:'%.4f'%x)
    def add_Data(self,Data,Tag='Sample'):
        self.DataBase.at[len(self.DataBase),'Data']=Data
        self.DataBase.at[len(self.DataBase)-1,'Data_Name'] = Data.DataName
        self.DataBase.at[len(self.DataBase)-1,'Tag'] = Tag
        self.DataBase.sort_values(by='Tag',ascending=False,inplace=True)
        self.DataBase.reset_index(drop=True,inplace=True)
    def show_Data(self):
        print(self.DataBase.iloc[:,[1,2]])
    def get_param(self):
        print(self.AlignmentParam)
    def set_param(self,Name,Value):
        self.AlignmentParam[Name] = Value
    def RenewRefList(self):
        self.RefList=pd.DataFrame(columns=['m/z','Int','RT','MS_List','RT_List','MS2_Int','MS2_MZ','SimilarityScore','SC_List'])
        self.RefList['m/z'] = self.RefList['m/z'].map(lambda x:'%.4f'%x)
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
        def RI_to_RT(RI,RIIS):
            n_place = np.where(RIIS['C']<=RI/100)[0]
            n1_place = np.where(RIIS['C']>RI/100)[0]
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
            RT = (RI/100-C_number)*(t_n1-t_n)+t_n
            return RT
        for i in range(len(self.DataBase)):
            self.RefList[self.DataBase.at[i,'Data_Name']] = 0
        for i in range(len(self.DataBase)):
            if self.DataBase.at[i,'Tag'] != 'Blank':
                self.DataBase.at[i,'Data'].Final_Peak_Detect.sort_values('AverageMZ',ascending=(False),inplace=True)
                self.DataBase.at[i,'Data'].Final_Peak_Detect.reset_index(drop=True,inplace=True)
                self.RefList.sort_values('m/z',ascending=(True),inplace=True)
                self.RefList.reset_index(drop=True,inplace=True)
                temp_Data = self.DataBase.at[i,'Data'].Final_Peak_Detect
                init_RefMZ = pd.DataFrame(columns=['m/z','MZList','RowIndex'])
                init_SampleMZ = pd.DataFrame(columns=['m/z','MZList','RowIndex'])
                if self.AlignmentParam['RI_Alignment'] == False:
                    Sample_Time = temp_Data.loc[:,'RT']
                else:
                    Sample_Time = temp_Data.loc[:,'RI']
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
                self.Label_process_sub.setText('Align '+str(i+1)+' / '+str(len(self.DataBase)))
                self.process_bar.setMaximum(len(init_SampleMZ))
                self.process_bar.setValue(0)
                QApplication.processEvents()
                add_MZ = []
                add_Int = []
                add_RT = []
                add_MS_List = []
                add_RT_List = []
                add_SimilarityScore = []
                add_SC_List = []
                for ii in range(len(init_SampleMZ)):
                    self.process_bar.setValue(self.process_bar.value()+1)
                    QApplication.processEvents()
                    Sample_MZ = init_SampleMZ.at[ii,'m/z']
                    match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
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
                            add_RT.append(Sample_Time[i_add])
                            add_MS_List.append([temp_Data.at[i_add,'AverageMZ']])
                            add_RT_List.append([Sample_Time[i_add]])
                            add_SimilarityScore.append(temp_Data.at[i_add,'SimilarityScore'])
                            add_SC_List.append([temp_Data.at[i_add,'SimilarityScore']])
                        continue
                    Score_matrix= np.zeros([len(Ref_Index),len(init_SampleMZ.at[ii,'RowIndex'])])
                    for i_Ref in range(len(Ref_Index)):
                        for i_Sample in range(len(init_SampleMZ.at[ii,'RowIndex'])):
                            i_Ref_MZ = self.RefList.at[Ref_Index[i_Ref],'m/z']
                            i_Sample_MZ = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sample],'AverageMZ']
                            i_Ref_Time = self.RefList.at[Ref_Index[i_Ref],'RT']
                            i_Sample_Time = Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sample]]
                            if abs(i_Ref_MZ-i_Sample_MZ)/i_Ref_MZ<self.AlignmentParam['MZ_Tor'] and abs(i_Ref_Time-i_Sample_Time)/i_Ref_Time<self.AlignmentParam['RT_Tor']:
                                Score_matrix[i_Ref,i_Sample] = self.AlignmentParam['A']*np.exp(-0.5*((i_Sample_Time-i_Ref_Time)/(i_Ref_Time*self.AlignmentParam['RT_Tor']))**2)+(1-self.AlignmentParam['A'])*np.exp(-0.5*((i_Sample_MZ-i_Ref_MZ)/(i_Sample_MZ*self.AlignmentParam['MZ_Tor']))**2)
                    Sm_row,Sm_col = linear_sum_assignment(Score_matrix,True)
                    for i_Sm_row,i_Sm_col in zip(Sm_row,Sm_col):
                        if Score_matrix[i_Sm_row,i_Sm_col] != 0:
                            self.RefList.at[Ref_Index[i_Sm_row],self.DataBase.at[i,'Data_Name']] = temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'Int']
                            self.RefList.at[Ref_Index[i_Sm_row],'MS_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                            self.RefList.at[Ref_Index[i_Sm_row],'RT_List'].append(Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                            self.RefList.at[Ref_Index[i_Sm_row],'SC_List'].append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'SimilarityScore'])
                            self.RefList.at[Ref_Index[i_Sm_row],'RT'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'RT_List'])
                            self.RefList.at[Ref_Index[i_Sm_row],'m/z'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'MS_List'])
                            self.RefList.at[Ref_Index[i_Sm_row],'SimilarityScore'] = np.mean(self.RefList.at[Ref_Index[i_Sm_row],'SC_List'])
                        else:
                            add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ'])
                            add_Int.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'Int'])
                            add_RT.append(Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]])
                            add_SimilarityScore.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'SimilarityScore'])
                            add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'AverageMZ']])
                            add_RT_List.append([Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col]]])
                            add_SC_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_Sm_col],'SimilarityScore']])
                    miss_col = list(filter(lambda x:x not in Sm_col,range(len(init_SampleMZ.at[ii,'RowIndex']))))
                    for i_miss in miss_col:
                        add_MZ.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ'])
                        add_Int.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'Int'])
                        add_RT.append(Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_miss]])
                        add_SimilarityScore.append(temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'SimilarityScore'])
                        add_MS_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'AverageMZ']])
                        add_RT_List.append([Sample_Time[init_SampleMZ.at[ii,'RowIndex'][i_miss]]])
                        add_SC_List.append([temp_Data.at[init_SampleMZ.at[ii,'RowIndex'][i_miss],'SimilarityScore']])
                add_RefList = pd.DataFrame({'m/z': add_MZ,'RT': add_RT,self.DataBase.at[i, 'Data_Name']: add_Int,
                                            'MS_List':add_MS_List,'RT_List':add_RT_List,'SimilarityScore':add_SimilarityScore,
                                            'SC_List':add_SC_List})
                self.RefList = pd.concat([self.RefList, add_RefList])
                self.RefList.reset_index(drop=True, inplace=True)
                self.RefList = self.RefList.fillna(0)
        Sample_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Sample'])
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        self.RefList['max_SampleInt'] = self.RefList.apply(lambda x:max(x[Sample_Name]),axis=1)
        self.RefList['Int'] = self.RefList.apply(lambda x:np.mean(x[Sample_Name]),axis=1)
        if self.AlignmentParam['RI_Alignment'] == True:
            self.RefList['RI'] = self.RefList['RT'].apply(lambda x:x)
            self.RefList['RT'] = self.RefList['RT'].apply(lambda x:RI_to_RT(x,self.DataBase['Data'][0].RIIS))
        if len(Blank_Name) > 0:
            for i in Blank_Name:
                temp_Blank_Data = self.DataBase[self.DataBase['Data_Name']==i].copy()
                temp_Blank_Data.reset_index(drop=True,inplace=True)
                temp_Blank_Data = temp_Blank_Data['Data'][0]
                for ii in range(len(self.RefList)):
                    [Blank_RT_List,Blank_Int_List] = temp_Blank_Data.ExtractDataPoint(self.RefList.at[ii,'m/z'],self.RefList.at[ii,'RT']+3,self.RefList.at[ii,'RT']-3,smooth_index=0)
                    self.RefList[i][ii] = max(1,max(Blank_Int_List))
            self.RefList['max_BlankInt'] = self.RefList.apply(lambda x:np.mean(x[Blank_Name]) if max(x[Blank_Name])>0 else 1,axis=1)
    def BlankFilter(self):
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        if len(Blank_Name) > 0:
            self.RefList.drop(self.RefList[self.RefList['max_SampleInt']/self.RefList['max_BlankInt']<self.AlignmentParam['Threshold']].index,inplace=True)
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
        