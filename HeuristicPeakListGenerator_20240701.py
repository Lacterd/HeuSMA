#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:30:30 2024

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

class HeuristicProducerUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
        # 设置标题
        self.resize(600,700)
        self.centerWindow()
        
        # 设置参数字典
        self.Data_params={}
        self.Calculate_params={'MS1_Tor':0.000010,'smooth':5,'min_Int':9000,'Points':17}
        
    def centerWindow(self):
        screen = QDesktopWidget().screenGeometry()
        size  = self.geometry()
        LeftValue  = int((screen.width()-size.width())/2)
        TopValue = int((screen.height()-size.height())/2)
        self.move(LeftValue,TopValue)
        
    def initUI(self):
        globallayout = QVBoxLayout()
        # Data import
        Data_import_Widget = QWidget()
        Data_import_Layout = QGridLayout()
        self.Label_SampleSelect = QLabel('Select data')
        self.Label_RIISSelect = QLabel('Select Calibrants data')
        self.Label_Params_Select = QLabel('Set params')
        self.TextBrowser_SampleSelect = QTableWidget()
        self.TextBrowser_SampleSelect.setColumnCount(2)
        self.TextBrowser_SampleSelect.setHorizontalHeaderLabels(['Data','Peak List'])
        self.TextBrowser_SampleSelect.setColumnWidth(0, 230)
        self.TextBrowser_SampleSelect.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.TextBrowser_SampleSelect.horizontalHeader().setStretchLastSection(True)
        self.TextBrowser_RIISSelect = QLabel('*')
        self.PushButton_SampleSelect = QPushButton('Select')
        self.PushButton_RIISSelect = QPushButton('Select')
        self.PushButton_SampleSelect.setToolTip('Select sample *.mzML documents')
        self.PushButton_RIISSelect.setToolTip('Select RIIS *.xlsx or *.xls document')
        # Slot
        self.PushButton_SampleSelect.clicked.connect(self.SelectSampleFile)
        self.PushButton_RIISSelect.clicked.connect(self.SelectRIISFile)
        Data_import_Layout.addWidget(self.Label_SampleSelect,1,0)
        Data_import_Layout.addWidget(self.TextBrowser_SampleSelect, 2, 0) # row 1， column 0
        Data_import_Layout.addWidget(self.PushButton_SampleSelect, 2, 1) # row 1， column 1
        Data_import_Layout.addWidget(self.Label_RIISSelect,5,0)
        Data_import_Layout.addWidget(self.TextBrowser_RIISSelect, 6, 0)
        Data_import_Layout.addWidget(self.PushButton_RIISSelect, 6, 1) # 行，列，行高，列宽
        Data_import_Widget.setLayout(Data_import_Layout)
        Params_setting_Widget = QWidget()
        Params_setting_Layout = QGridLayout()
        # Union
        self.Lable_MS_Tor = QLabel('MS1 Tolerance(ppm)')
        self.Lable_RT_Tor = QLabel('RI Tolerance(%)')
        self.LineEdit_MS_Tor = QLineEdit('10')
        self.LineEdit_RT_Tor = QLineEdit('2')
        self.LineEdit_MS_Tor.setToolTip('Under this threshold will be recognized as same m/z')
        self.LineEdit_RT_Tor.setToolTip('Under this threshold will be recognized as same RT')
        Params_setting_Layout.addWidget(self.Label_Params_Select,0,0)
        Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.LineEdit_MS_Tor,1,1)
        Params_setting_Layout.addWidget(self.Lable_RT_Tor,1,2)
        Params_setting_Layout.addWidget(self.LineEdit_RT_Tor,1,3)
        Params_setting_Widget.setLayout(Params_setting_Layout)
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
        self.setWindowTitle('Heuristic Peak List Generator')

    def SelectSampleFile(self):
        FileName,FileType = QFileDialog.getOpenFileNames(self,"选取文件",os.getcwd(),"mzML Files(*.mzML)")
        for i in FileName:
            name_begin = i.rfind('/')
            if name_begin == -1:
                name_begin = FileName.rfind('\\')
            name_end = len(i)
            LineEdit_Data= QLineEdit(i[name_begin+1:name_end])
            LineEdit_Data.textChanged.connect(self.DataChange)
            LineEdit_PeakList= QLineEdit(i[name_begin+1:name_end-5]+'.xlsx')
            LineEdit_PeakList.textChanged.connect(self.PeakListChange)
            self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,0,LineEdit_Data)
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,LineEdit_PeakList)
            self.Data_params[i[name_begin+1:name_end-5]] = {'Name':i[name_begin+1:name_end-5],'Path':i,'PeakListPath':i[0:name_end-5]+'.xlsx','LineEdit_Data':LineEdit_Data,'LineEdit_PeakList':LineEdit_PeakList}
            for i in self.Data_params.keys():
                if sys.platform == 'darwin':
                    name_begin = self.Data_params[i]['Path'].rfind('/')
                    self.filepath_title = self.Data_params[i]['Path'][0:name_begin]+'/'
                elif sys.platform == 'win':
                    name_begin = self.Data_params[i]['Path'].rfind('\\')
                    self.filepath_title = self.Data_params[i]['Path'][0:name_begin]+'\\'
                break

    def SelectRIISFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"选取文件",os.getcwd(),"Excel Files(*.xlsx)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_RIISSelect.setText(FileName[name_begin+1:name_end])
        self.RIISPath = FileName
    
    def processbar_fresh(self,*arg):
        self.process_bar.setValue(self.process_bar.value()+1)
        QApplication.processEvents()
        
    def DataChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['LineEdit_Data']:
               self.Data_params[i]['Path'] = self.filepath_title + self.Data_params[i]['LineEdit_Data'].text()
    
    def PeakListChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['LineEdit_PeakList']:
               self.Data_params[i]['PeakListPath'] = self.filepath_title + self.Data_params[i]['LineEdit_PeakList'].text()
        
    def Run(self):
        self.SampleData = {}
        for i in self.Data_params.keys():
            self.SampleData[i] = self.Data_params[i]
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
            self.SampleData[i]['Data'].set_RIIS(self.RIISPath)
            self.SampleData[i]['Data'].Final_Peak_Detect = pd.read_excel(self.SampleData[i]['PeakListPath'])
            self.SampleData[i]['Data'].Final_Peak_Detect.rename(columns={'m/z':'AverageMZ'},inplace=True)
            self.SampleData[i]['Data'].Calculate_RI()
            self.SampleData[i]['Data'].Final_Peak_Detect['RI_left'] = self.SampleData[i]['Data'].Final_Peak_Detect['RI'].apply(lambda x:x*0.98)
            self.SampleData[i]['Data'].Final_Peak_Detect['RI_right'] = self.SampleData[i]['Data'].Final_Peak_Detect['RI'].apply(lambda x:x*1.02)
            self.SampleData[i]['Data'].Final_Peak_Detect['MZ_List'] = self.SampleData[i]['Data'].Final_Peak_Detect['AverageMZ'].apply(lambda x:[x])
            self.SampleData[i]['Data'].Final_Peak_Detect['RI_List'] = self.SampleData[i]['Data'].Final_Peak_Detect['RI'].apply(lambda x:[x])
        Alltemp_HPL=pd.DataFrame(columns=['AverageMZ','RT','MZ_List','RT_List'])
        Alltemp_HPL['AverageMZ'] = Alltemp_HPL['AverageMZ'].map(lambda x:'%.4f'%x)
        sub_HPL_dict = {}
        for i in self.SampleData.keys():
            sub_HPL_dict[i] = self.SampleData[i]['Data'].Final_Peak_Detect.copy()
        MS_Tor = float(self.LineEdit_MS_Tor.text())/1000000
        RT_Tor = float(self.LineEdit_RT_Tor.text())/100
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
                    if abs(init_RefMZ.at[len(init_RefMZ)-1,'AverageMZ']-Alltemp_HPL.at[ii,'AverageMZ'])/Alltemp_HPL.at[ii,'AverageMZ']<MS_Tor:
                        init_RefMZ.at[len(init_RefMZ)-1,'MZList'].append(Alltemp_HPL.at[ii,'AverageMZ'])
                        init_RefMZ.at[len(init_RefMZ)-1,'RowIndex'].append(ii)
                        init_RefMZ.at[len(init_RefMZ)-1,'AverageMZ']=np.mean(init_RefMZ.at[len(init_RefMZ)-1,'MZList'])
                    else:
                        init_RefMZ.loc[len(init_RefMZ)] = [Alltemp_HPL.at[ii,'AverageMZ'],[Alltemp_HPL.at[ii,'AverageMZ']],[ii]]
            for ii in range(len(temp_Data)):
                if len(init_SampleMZ)==0:
                    init_SampleMZ.loc[len(init_SampleMZ)] = [temp_Data.at[ii,'AverageMZ'],[temp_Data.at[ii,'AverageMZ']],[ii]]
                else:
                    if abs(init_SampleMZ.at[len(init_SampleMZ)-1,'AverageMZ']-temp_Data.at[ii,'AverageMZ'])/temp_Data.at[ii,'AverageMZ']<MS_Tor:
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
                match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'AverageMZ']-Sample_MZ)/Sample_MZ<MS_Tor,range(len(init_RefMZ))))
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
                        if abs(i_Ref_MZ-i_Sample_MZ)/i_Ref_MZ<MS_Tor and abs(i_Ref_Time-i_Sample_Time)/i_Ref_Time<RT_Tor:
                            Score_matrix[i_Ref,i_Sample] = 0.5*np.exp(-0.5*((i_Sample_Time-i_Ref_Time)/(i_Ref_Time*RT_Tor))**2)+(1-0.5)*np.exp(-0.5*((i_Sample_MZ-i_Ref_MZ)/(i_Sample_MZ*MS_Tor))**2)
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
        Alltemp_HPL['RI_left'] = Alltemp_HPL['RI'].apply(lambda x:x*(1-RT_Tor))
        Alltemp_HPL['RI_right'] = Alltemp_HPL['RI'].apply(lambda x:x*(1+RT_Tor))
        Alltemp_HPL.sort_values('AverageMZ',ascending=(True),inplace=True)
        Alltemp_HPL.reset_index(drop=True,inplace=True)
        self.Heuristic_normal_list = Alltemp_HPL.copy()
        if sys.platform == 'darwin':
            name_begin = self.SampleData[list(self.SampleData.keys())[0]]['Path'].rfind('/')
            file_output = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'/HeuristicList'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx'
            pickle_path = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'/HeuristicList'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl'
        elif sys.platform == 'win':
            name_begin = self.SampleData[list(self.SampleData.keys())[0]]['Path'].rfind('\\')
            file_output = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'\\HeuristicList'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx'
            pickle_path = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'\\HeuristicList'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl'
        self.Heuristic_normal_list.rename(columns={'AverageMZ':'Average m/z'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI':'Average RI'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'MZ_List':'m/z value in each sample'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI_List':'RI value in each sample'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI_left':'RI left boundary'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI_right':'RI right boundary'},inplace=True)
        self.Heuristic_normal_list.to_excel(file_output,index=False)
        with open(pickle_path,'wb') as f:
            pickle.dump(self.Heuristic_normal_list,f)
        self.Label_process_sub.setText('Finished')

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
        '''
        if plt_for_test == True:
            plt.figure()
            plt.title(str(MZ))
            plt.plot(RT_List,Int_List)
            plt.figure()
        '''
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


if __name__ == '__main__':
    mp.freeze_support()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('./HeuristicList.ico'))
    #main = FirstMainWindow()
    main_HPLG = HeuristicProducerUI()
    main_HPLG.show()
    sys.exit(app.exec_())
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        