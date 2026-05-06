#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:30:30 2024

@author: lacter
"""

import pandas as pd
import time
from progress.bar import Bar
from pyopenms import MSExperiment,MzMLFile,MzXMLFile
import glob
import math
import numpy as np
import bisect
import re
import os
from scipy.optimize import linear_sum_assignment
from scipy.spatial import KDTree
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
import pyopenms

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
        self.TextBrowser_SampleSelect.setColumnCount(3) #修改
        self.TextBrowser_SampleSelect.setHorizontalHeaderLabels(['Data','Gradient','Blank']) #修改
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
        self.Lable_RI_Tor = QLabel('RI Tolerance(%)')
        self.Lable_Int_Tor = QLabel('Int Tolerance(%)')
        self.Lable_Parallel_Tor = QLabel('Detection frequency')
        self.Lable_RT_Tor = QLabel('RT Tolerance(min)')
        self.Lable_Point = QLabel('Min scan points')
        self.Lable_SN = QLabel('Signal/Noise')
        self.Lable_SB = QLabel('Signal/Blank')
        self.LineEdit_MS_Tor = QLineEdit('10')
        self.LineEdit_RI_Tor = QLineEdit('2')
        self.LineEdit_Int_Tor = QLineEdit('10000')
        self.LineEdit_Parallel_Tor = QLineEdit('0.5')
        self.LineEdit_Point = QLineEdit('17')
        self.LineEdit_SN = QLineEdit('10')
        self.LineEdit_SB = QLineEdit('3')
        self.LineEdit_RT_Tor = QLineEdit('0.1')
        self.LineEdit_MS_Tor.setToolTip('Under this threshold will be recognized as same m/z')
        self.LineEdit_RI_Tor.setToolTip('Under this threshold will be recognized as same RI')
        self.LineEdit_Int_Tor.setToolTip('Minimize intensity')
        self.LineEdit_Point.setToolTip('Scan points under this threshold will be deleted')
        self.LineEdit_SN.setToolTip('Under this threshold will be deleted')
        self.LineEdit_SB.setToolTip('Under this threshold will be deleted')
        self.LineEdit_Parallel_Tor.setToolTip('Peaks should appear in samples beyond this ratio')
        self.LineEdit_RT_Tor.setToolTip('Under this threshold will be recognized as same RT')
        Params_setting_Layout.addWidget(self.Label_Params_Select,0,0)
        Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.LineEdit_MS_Tor,1,1)
        Params_setting_Layout.addWidget(self.Lable_RI_Tor,1,2)
        Params_setting_Layout.addWidget(self.LineEdit_RI_Tor,1,3)
        Params_setting_Layout.addWidget(self.Lable_Int_Tor,2,0)
        Params_setting_Layout.addWidget(self.LineEdit_Int_Tor,2,1)
        Params_setting_Layout.addWidget(self.Lable_Parallel_Tor,2,2)
        Params_setting_Layout.addWidget(self.LineEdit_Parallel_Tor,2,3)
        Params_setting_Layout.addWidget(self.Lable_RT_Tor,3,0)
        Params_setting_Layout.addWidget(self.LineEdit_RT_Tor,3,1)
        Params_setting_Layout.addWidget(self.Lable_Point,3,2)
        Params_setting_Layout.addWidget(self.LineEdit_Point,3,3)
        Params_setting_Layout.addWidget(self.Lable_SB,4,0)
        Params_setting_Layout.addWidget(self.LineEdit_SB,4,1)
        Params_setting_Layout.addWidget(self.Lable_SN,4,2)
        Params_setting_Layout.addWidget(self.LineEdit_SN,4,3)
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
            LineEdit_Gradient= QLineEdit('') 
            LineEdit_Gradient.textChanged.connect(self.GradientChange)
            LineEdit_PeakBlank= QLineEdit('') 
            LineEdit_PeakBlank.textChanged.connect(self.PeakBlank)
            self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,0,LineEdit_Data)
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,LineEdit_Gradient) #修改
            self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,2,LineEdit_PeakBlank) #修改
            self.Data_params[i[name_begin+1:name_end-5]] = {'Name':i[name_begin+1:name_end-5],'Path':i,'LineEdit_Gradient':LineEdit_Gradient,'LineEdit_Data':LineEdit_Data,'LineEdit_PeakBlank':LineEdit_PeakBlank,'Gradient':'','PeakBlank':''}
        for i in self.Data_params.keys():
            name_begin = self.Data_params[i]['Path'].rfind('/')
            if name_begin > 0:
                self.filepath_title = self.Data_params[i]['Path'][0:name_begin]+'/'
            else:
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
               input_text = self.Data_params[i]['LineEdit_Data'].text()
               self.Data_params[i]['Path'] = self.filepath_title + input_text
    
    def GradientChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['LineEdit_Gradient']:
               input_text = self.Data_params[i]['LineEdit_Gradient'].text()
               self.Data_params[i]['Gradient'] = 'Gradient-'+ input_text
    
    def PeakBlank(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['LineEdit_PeakBlank']:
               self.Data_params[i]['PeakBlank'] = self.filepath_title + self.Data_params[i]['LineEdit_PeakBlank'].text()
        
    def Run(self):
        for i in self.Data_params.keys():
            input_text = self.Data_params[i]['LineEdit_Gradient'].text()
            self.Data_params[i]['Gradient'] = input_text
        for i in self.Data_params.keys():
            if len(self.Data_params[i]['LineEdit_PeakBlank'].text()) > 0:
                self.Data_params[i]['PeakBlank'] = self.filepath_title + self.Data_params[i]['LineEdit_PeakBlank'].text()
            else:
                self.Data_params[i]['PeakBlank'] = ''
        self.SampleData = {}
        self.GradientList = []
        for i in self.Data_params.keys():
            self.SampleData[i] = self.Data_params[i]
            self.GradientList.append(self.Data_params[i]['Gradient'])
        self.GradientList = list(set(self.GradientList))
        self.Label_process_sub.setText('Load Data')
        self.process_bar.setMaximum(len(self.SampleData.keys()))
        self.process_bar.setValue(0)
        ''' -Sample- '''
        self.Label_process_sub.setText('Load Data')
        self.process_bar.setMaximum(len(self.SampleData.keys()))
        self.process_bar.setValue(0)
        pool = mp.Pool(os.cpu_count()-2)
        self.pool_result = {}
        QApplication.processEvents()
        for i in self.SampleData.keys():
            self.pool_result[i]=pool.apply_async(pool_HPLG,args=(self.SampleData[i]['Path'],float(self.LineEdit_MS_Tor.text())/1000000,float(self.LineEdit_RI_Tor.text())/100,float(self.LineEdit_RT_Tor.text())*60,int(self.LineEdit_Int_Tor.text()),int(self.LineEdit_Point.text()),self.RIISPath,int(self.LineEdit_SN.text()),self.SampleData[i]['PeakBlank'],int(self.LineEdit_SB.text()),self.Data_params[i]['Gradient'],))
        pool.close()
        pool.join()
        for i in self.SampleData.keys():
            self.SampleData[i]['Final_Peak_Detect'] = self.pool_result[i].get().copy()
        temp_key = list(self.SampleData.keys())[0]
        temp_EMZDP = EazyMZDataProcess(self.SampleData[temp_key]['Path'])
        self.Polarity = temp_EMZDP.get_param('Polarity')
        ''' -- '''
        name_begin = self.SampleData[list(self.SampleData.keys())[0]]['Path'].rfind('/')
        if name_begin > 0:
            self.Path_title = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin+1]
            file_output = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'/HeuristicList-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx'
            pickle_path = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'/HeuristicList-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl'
        else:
            name_begin = self.SampleData[list(self.SampleData.keys())[0]]['Path'].rfind('\\')
            self.Path_title = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin+2]
            file_output = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'\\HeuristicList-'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx'
            pickle_path = self.SampleData[list(self.SampleData.keys())[0]]['Path'][0:name_begin]+'\\HeuristicLis-t'+self.Polarity+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.pkl'
        ''' in gradient alignment '''
        Align_Graident = {}
        for i in self.GradientList:
            temp_Align = DataAlignment()
            temp_Align.AlignmentParam['RI_Alignment'] = True
            temp_Align.AlignmentParam['Miss_Filter'] = float(self.LineEdit_Parallel_Tor.text())
            temp_Align.AlignmentParam['RI_Tor'] = float(self.LineEdit_RI_Tor.text())/100
            for ii in self.SampleData.keys():
                if self.SampleData[ii]['Gradient'] == i:
                    temp_Align.add_Data_new(self.SampleData[ii],Tag='Sample')
            temp_Align.RenewRefList_new()
            temp_Align.Filter_MissingValue(WhetherDel=True)
            temp_Align.RefList.rename(columns={'m/z':'AverageMZ'},inplace=True)
            Align_Graident[i] = {}
            Align_Graident[i]['Final_Peak_Detect'] = temp_Align.RefList.copy()
            Align_Graident[i]['Name'] = i
            temp_Align.RefList.to_excel(self.Path_title+self.Polarity+'-'+i+'-in-gradient-alignment'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx')
        temp_Align = DataAlignment()
        temp_Align.AlignmentParam['RI_Alignment'] = True
        temp_Align.AlignmentParam['RI_Tor'] = float(self.LineEdit_RI_Tor.text())/100
        for i in self.GradientList:
            temp_Align.add_Data_new(Align_Graident[i],Tag='Sample')
        temp_Align.RenewRefList_new()
        temp_Align.RefList.rename(columns={'m/z':'AverageMZ'},inplace=True)
        temp_Align.RefList['RI_left'] = temp_Align.RefList['RI'].apply(lambda x:x*(1-float(self.LineEdit_RI_Tor.text())/100))
        temp_Align.RefList['RI_right'] = temp_Align.RefList['RI'].apply(lambda x:x*(1+float(self.LineEdit_RI_Tor.text())/100))
        temp_Align.RefList.sort_values('AverageMZ',ascending=(True),inplace=True)
        temp_Align.RefList.reset_index(drop=True,inplace=True)
        temp_Align.RefList.drop('RT_List',axis=1,inplace=True)
        temp_Align.RefList.drop('MS2_Int',axis=1,inplace=True)
        temp_Align.RefList.drop('MS2_MZ',axis=1,inplace=True)
        temp_Align.RefList.drop('max_SampleInt',axis=1,inplace=True)
        temp_Align.RefList.drop('max_BlankInt',axis=1,inplace=True)
        temp_Align.RefList.drop('Int',axis=1,inplace=True)
        self.Heuristic_normal_list = temp_Align.RefList.copy()
        with open(pickle_path,'wb') as f:
            pickle.dump(self.Heuristic_normal_list,f)
        self.Heuristic_normal_list.rename(columns={'AverageMZ':'Average m/z'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI':'Average RI'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'MS_List':'m/z value in each sample'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI_List':'RI value in each sample'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI_left':'RI left boundary'},inplace=True)
        self.Heuristic_normal_list.rename(columns={'RI_right':'RI right boundary'},inplace=True)
        self.Heuristic_normal_list.to_excel(file_output,index=False)
        self.Label_process_sub.setText('Finished')

def pool_HPLG(Path,MS1_Tor,RI_Tor,RT_Tor,min_Int,Points,RIISPath,SN,BlankPath,SB,Gradient_Name):
    temp_EMZDP = EazyMZDataProcess(Path)
    temp_EMZDP.set_param('MS1_Tor',MS1_Tor)
    temp_EMZDP.set_param('RI_Tor',RI_Tor)
    temp_EMZDP.set_param('RT_Tor',RT_Tor)
    temp_EMZDP.set_param('min_Int',min_Int)
    temp_EMZDP.set_param('Points',Points)
    RIIS_temp = pd.read_excel(RIISPath)
    if Gradient_Name in list(RIIS_temp.keys()):
        RIISPath = pd.DataFrame({'C':RIIS_temp.loc[:,'C'],'RT':RIIS_temp.loc[:,Gradient_Name]})
    temp_EMZDP.set_RIIS(RIISPath)
    File_path = glob.glob(Path.replace('.mzML','.pkl'))
    if len(File_path) > 0:
        temp_EMZDP.load_FPD(File_path[0])
    else:
        temp_EMZDP.Final_Peak_Detect = temp_EMZDP.detect_Peak()
        temp_EMZDP.Calculate_RI()
        temp_EMZDP.save_FPD(Path.replace('.mzML','-Origin.pkl'))
        if len(BlankPath)>0:
            temp_EMZDP.add_Blank(BlankPath)
            temp_EMZDP.Calculate_SB(drop=True,limit=SB)
        temp_EMZDP.Calculate_SN(drop=True,Threshold=SN)
        temp_EMZDP.Deconvolution(WhetherDel=True)
        temp_EMZDP.save_FPD(Path.replace('.mzML','.pkl'))
    return temp_EMZDP.Final_Peak_Detect.copy()

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
        self.__param = {'MS1_Tor':0.000010,'RT_Tor':6,'min_Int':10000,'min_RT':60,'min_RT_width':6,'max_Noise':2000,
                        'RI_Tor':0.02,'Deconvolution':False,'FeatureDetectPlot':3,'MergeRule':'Intersection',
                        'UpDown_gap':10,'saveAutoList':False,'smooth':5,'Points':17,'DeconvolutionSimilarityScore':0.98,
                        'assign MS2':True,'min_MZ':100}
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
                #temp_range = list(filter(lambda x:MZ_temp[x]<=Pre_temp,range(len(MZ_temp)))) # and Int_temp[x]>1000
                '''
                temp_range = list(filter(lambda x:len(np.where(abs(MZ_temp[temp_range]-MZ_temp[x])/MZ_temp[x]<self.__param['MS1_Tor']))==1 or
                                         Int_temp[x]==max(Int_temp[np.where(abs(MZ_temp[temp_range]-MZ_temp[x])/MZ_temp[x]<self.__param['MS1_Tor'])]),temp_range))
                '''
                if True:#len(temp_range)>0
                    #MZ_temp = MZ_temp[temp_range]
                    #Int_temp = Int_temp[temp_range]
                    self.MS2_Pre.append(Pre_temp)
                    self.MS2_RT_List.append(i.getRT())
                    #self.MS2_MZ_List.append(MZ_temp[temp_range])
                    #self.MS2_Int_List.append(Int_temp[temp_range])
                    self.MS2_MZ_List.append(MZ_temp)
                    self.MS2_Int_List.append(Int_temp)
                    #self.MS2_RelInt.append(max(Int_temp[temp_range]))
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
        self.MS2_Data = {'Scan_Time':self.MS2_RT_List,'Pre_MZ':self.MS2_Pre,'MZ_List':self.MS2_MZ_List,'Int_List':self.MS2_Int_List}
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
        self.Final_Peak_Detect = pd.DataFrame(columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])

    def add_0(x):
        x.append(0)
        return x
    def add_Blank(self,DataPath):
        self.BlankData = MSExperiment()
        if DataPath.endswith('mzML'):
            MzMLFile().load(DataPath,self.BlankData)
        elif DataPath.endswith('mzXML'):
            MzXMLFile().load(DataPath,self.BlankData)
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
    def assign_MS2(self,RI_swich=False):
        for i in self.MS2_Data.keys():  
            V_List_Num = list(self.MS2_Data['Pre_MZ'][int(len(self.MS2_Data['Pre_MZ'])/2):int(len(self.MS2_Data['Pre_MZ'])/2)+50])
            V_List_Count = []
            Count = 1
            for ii in range(1,len(V_List_Num)):
                if abs(V_List_Num[ii-1]-V_List_Num[ii])/V_List_Num[ii] <= self.get_param('MS1_Tor'):
                    Count = Count + 1
                else:
                    V_List_Count.append(Count)
                    Count = 1
            V_List_Count = pd.Series(V_List_Count[2:len(V_List_Count)-2])
            V_List = []
            for ii in range(V_List_Count.mode()[0]):
                V_List.append('CE'+str(ii+1))
            break
        self.MS2_Data['CE'] = V_List[0]
        for i in range(1,len(self.MS2_Data)):
            if self.MS2_Data['Pre_MZ'][i] == self.MS2_Data['Pre_MZ'][i-1]:
                V_Place = list(filter(lambda x:V_List[x]==self.MS2_Data['CE'][i-1],range(len(V_List))))[0]
                if V_Place == len(V_List)-1:
                    self.MS2_Data.at[i,'CE'] = V_List[0]
                else:
                    self.MS2_Data.at[i,'CE'] = V_List[V_Place+1]
            else:
                self.MS2_Data.at[i,'CE'] = V_List[0]
        if 'MS2_MZ' not in self.Final_Peak_Detect.keys():
            self.Final_Peak_Detect['MS2_MZ'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:[])
            self.Final_Peak_Detect['MS2_Int'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:[])
        self.Final_Peak_Detect['MS2_MZ'] = self.Final_Peak_Detect['MS2_MZ'].astype('object')
        self.Final_Peak_Detect['MS2_Int'] = self.Final_Peak_Detect['MS2_Int'].astype('object')
        if RI_swich == False:
            bar = Bar('Assign MS2',max=len(self.Final_Peak_Detect))
            data_points = np.column_stack((self.MS2_Data['Scan_Time'],self.MS2_Data['Pre_MZ']))
            radius = self.get_param('RT_Tor')
            kdtree = KDTree(data_points)
            for i,RT,MZ in zip(range(len(self.Final_Peak_Detect)),self.Final_Peak_Detect['RT'],self.Final_Peak_Detect['AverageMZ']):
                bar.next()
                Match_List = kdtree.query_ball_point([RT,MZ], r=radius)
                Match_List = list(filter(lambda x:abs(self.MS2_Data.at[x,'Pre_MZ']-MZ)/MZ<self.get_param('MS1_Tor'),Match_List))
                if len(Match_List)>0:
                    Match_List = list(filter(lambda x:abs(RT-self.MS2_Data.loc[x,'Scan_Time'])==min(abs(np.array(self.MS2_Data.loc[Match_List,'Scan_Time'])-RT)),Match_List))[0]
                    CE_Index = list(filter(lambda x:V_List[x]==self.MS2_Data.at[Match_List,'CE'],range(len(V_List))))[0]
                    if Match_List-CE_Index > 0 and Match_List+len(V_List)-CE_Index < len(self.MS2_Data):
                        Valified_Index_List = list(range(Match_List-CE_Index,Match_List+len(V_List)-CE_Index))
                    elif Match_List-CE_Index < 0:
                        Valified_Index_List = list(range(0,Match_List+len(V_List)-CE_Index))
                    elif Match_List+len(V_List)-CE_Index > len(self.MS2_Data):
                        Valified_Index_List = list(range(Match_List-CE_Index,len(self.MS2_Data)))
                    MZ_List = []
                    Int_List = []
                    for ii in Valified_Index_List:
                        if abs(self.MS2_Data.at[ii,'Pre_MZ'] - MZ) / MZ<self.get_param('MS1_Tor') and len(self.MS2_Data.at[ii,'Int_List'])>0:
                            if len(MZ_List) == 0:
                                MZ_List = list(self.MS2_Data.at[ii,'MZ_List'][(self.MS2_Data.at[ii,'Int_List']>max(self.MS2_Data.at[ii,'Int_List'])*0.05) & (self.MS2_Data.at[ii,'MZ_List'] <= MZ*(1+self.get_param('MS1_Tor')))])
                                Int_List = list(self.MS2_Data.at[ii,'Int_List'][(self.MS2_Data.at[ii,'Int_List']>max(self.MS2_Data.at[ii,'Int_List'])*0.05) & (self.MS2_Data.at[ii,'MZ_List'] <= MZ*(1+self.get_param('MS1_Tor')))])
                            else:
                                min_limit = max(self.MS2_Data.at[ii,'Int_List'])*0.05
                                for iii in range(len(self.MS2_Data.at[ii,'MZ_List'])):
                                    if self.MS2_Data.at[ii,'Int_List'][iii] >= min_limit and self.MS2_Data.at[ii,'MZ_List'][iii] <= MZ*(1+self.get_param('MS1_Tor')):
                                        MS2_Merge_Index = range(bisect.bisect_left(MZ_List,self.MS2_Data.at[ii,'MZ_List'][iii]*(1-self.get_param('MS1_Tor'))),bisect.bisect_right(MZ_List,self.MS2_Data.at[ii,'MZ_List'][iii]*(1+self.get_param('MS1_Tor'))))
                                        #MS2_Merge_Index = list(filter(lambda x:abs(MZ_List[x]-self.MS2_Data.at[ii,'MZ_List'][iii])/MZ_List[x]<self.get_param('MS1_Tor'),range(len(MZ_List))))
                                        if len(MS2_Merge_Index) == 0:
                                            MS2_Merge_Index = bisect.bisect_left(MZ_List,self.MS2_Data.at[ii,'MZ_List'][iii])
                                            MZ_List.insert(MS2_Merge_Index,self.MS2_Data.at[ii,'MZ_List'][iii])
                                            Int_List.insert(MS2_Merge_Index,self.MS2_Data.at[ii,'Int_List'][iii])
                                        elif Int_List[MS2_Merge_Index[0]] < self.MS2_Data.at[ii,'Int_List'][iii]:
                                            #MZ_List[MS2_Merge_Index[0]] = self.MS2_Data.at[ii,'MZ_List'][iii]
                                            Int_List[MS2_Merge_Index[0]] = self.MS2_Data.at[ii,'Int_List'][iii]
                    '''
                    temp_range = pd.DataFrame({'MZ_List':MZ_List,'Int_List':Int_List})
                    temp_range.sort_values(by='MZ_List',inplace=True)
                    temp_range.reset_index(drop=True,inplace=True)
                    self.Final_Peak_Detect.at[i,'MS2_MZ'] = list(temp_range['MZ_List'])
                    self.Final_Peak_Detect.at[i,'MS2_Int'] = list(temp_range['Int_List'])
                    '''
                    self.Final_Peak_Detect.at[i,'MS2_MZ'] = MZ_List
                    self.Final_Peak_Detect.at[i,'MS2_Int'] = Int_List
            bar.finish() 
        else:
            bar = Bar('Assign MS2',max=len(self.Final_Peak_Detect))
            data_points = np.column_stack((self.MS2_Data['Scan_RI'],self.MS2_Data['Pre_MZ']))
            radius = np.sqrt(self.get_param('MS1_Tor')**2 + self.get_param('RI_Tor')**2)
            kdtree = KDTree(data_points)
            for i,RI,MZ in zip(range(len(self.Final_Peak_Detect)),self.Final_Peak_Detect['RI'],self.Final_Peak_Detect['AverageMZ']):
                bar.next()
                Match_List = kdtree.query_ball_point([RI,MZ], r=radius)
                Match_List = list(filter(lambda x:abs(self.MS2_Data.at[x,'Pre_MZ']-MZ)/MZ<self.get_param('MS1_Tor'),Match_List))
                if len(Match_List)>0:
                    Match_List = list(filter(lambda x:abs(RI-self.MS2_Data.loc[x,'Scan_RI'])==min(abs(np.array(self.MS2_Data.loc[Match_List,'Scan_RI'])-RI)),Match_List))[0]
                    CE_Index = list(filter(lambda x:V_List[x]==self.MS2_Data.at[Match_List,'CE'],range(len(V_List))))[0]
                    if Match_List-CE_Index > 0 and Match_List+len(V_List)-CE_Index < len(self.MS2_Data):
                        Valified_Index_List = list(range(Match_List-CE_Index,Match_List+len(V_List)-CE_Index))
                    elif Match_List-CE_Index < 0:
                        Valified_Index_List = list(range(0,Match_List+len(V_List)-CE_Index))
                    elif Match_List+len(V_List)-CE_Index > len(self.MS2_Data):
                        Valified_Index_List = list(range(Match_List-CE_Index,len(self.MS2_Data)))
                    MZ_List = []
                    Int_List = []
                    for ii in Valified_Index_List:
                        if len(self.MS2_Data.at[ii,'Int_List'])>0: #abs(self.MS2_Data.at[ii,'Pre_MZ'] - MZ) / MZ<self.get_param('MS1_Tor')*2 and
                            if len(MZ_List) == 0:
                                MZ_List = list(self.MS2_Data.at[ii,'MZ_List'][self.MS2_Data.at[ii,'Int_List']>max(self.MS2_Data.at[ii,'Int_List'])*0.05])
                                Int_List = list(self.MS2_Data.at[ii,'Int_List'][self.MS2_Data.at[ii,'Int_List']>max(self.MS2_Data.at[ii,'Int_List'])*0.05])
                            else:
                                #min_limit = max(self.MS2_Data.at[ii,'Int_List'])*0.05
                                for iii in range(len(self.MS2_Data.at[ii,'MZ_List'])):
                                    #if self.MS2_Data.at[ii,'Int_List'][iii] >= min_limit:
                                        MS2_Merge_Index = list(filter(lambda x:abs(MZ_List[x]-self.MS2_Data.at[ii,'MZ_List'][iii])/MZ_List[x]<self.get_param('MS1_Tor'),range(len(MZ_List))))
                                        if len(MS2_Merge_Index) == 0:
                                            MZ_List.append(self.MS2_Data.at[ii,'MZ_List'][iii])
                                            Int_List.append(self.MS2_Data.at[ii,'Int_List'][iii])
                                        elif Int_List[MS2_Merge_Index[0]] < self.MS2_Data.at[ii,'Int_List'][iii]:
                                            MZ_List[MS2_Merge_Index[0]] = self.MS2_Data.at[ii,'MZ_List'][iii]
                                            Int_List[MS2_Merge_Index[0]] = self.MS2_Data.at[ii,'Int_List'][iii]
                    temp_range = pd.DataFrame({'MZ_List':MZ_List,'Int_List':Int_List})
                    temp_range.sort_values(by='MZ_List',inplace=True)
                    temp_range.reset_index(drop=True,inplace=True)
                    self.Final_Peak_Detect.at[i,'MS2_MZ'] = list(temp_range['MZ_List'])
                    self.Final_Peak_Detect.at[i,'MS2_Int'] = list(temp_range['Int_List'])
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
    def AdjustedCosineSimilarity(ExestList, NewList):
        if  len(ExestList) == 0 :
            return 0
        else:
            if len(ExestList) != len(NewList):
                print("CosineSimilarity error", len(ExestList), len(NewList))
                return 0
        ExestList = np.array(ExestList)-np.mean(ExestList)
        NewList = np.array(NewList)-np.mean(NewList)
        CosineSimilarityValue = ExestList.dot(NewList)/(np.linalg.norm(ExestList)*np.linalg.norm(NewList))
        return CosineSimilarityValue
        
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
        bar = Bar('Calculate_SB ',max=len(self.Final_Peak_Detect))
        for i in range(len(self.Final_Peak_Detect)):
            bar.next()
            RTL = self.Final_Peak_Detect.at[i,'RTList'][0]
            RTR = self.Final_Peak_Detect.at[i,'RTList'][-1]
            MZ = self.Final_Peak_Detect.at[i,'AverageMZ']
            [RT_List,Int_List] = self.ExtractBlankPoint(MZ,RTR,RTL,smooth_index=0)
            if len(Int_List)>0:
                self.Final_Peak_Detect.at[i,'BLANK'] = max(1,max(Int_List))
            else:
                self.Final_Peak_Detect.at[i,'BLANK'] = 1
        bar.finish()
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
        RT = (RI/100-C_number)*(t_n1-t_n)+t_n
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
        RI = 100*(C_number+(RT-t_n)/(t_n1-t_n))
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
            RI = 100*(C_number+(RT-t_n)/(t_n1-t_n))
            return RI
        gaussian_params = [1, 1.3, 1.6, 2]
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
            #Int_List = list(map(lambda x:sum(Int_List[x-2:x+2])/len(Int_List[x-2:x+2]) if x in range(2,len(Int_List)-2) else Int_List[x],range(len(Int_List))))
            Int_List = np.array(Int_List)
            Int_List[2:-2] = np.convolve(Int_List, np.ones(5)/5, mode='valid')
            Int_List = list(Int_List)
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
            #FD_P = list(filter(lambda x: FD_List[x] > FD_Median, range(1, len(FD_List))))
            #FD_N = list(filter(lambda x: FD_List[x] < FD_Median_N, range(1, len(FD_List))))
            FD_P = np.where(FD_List[1:] > FD_Median)[0] + 1
            FD_N = np.where(FD_List[1:] < FD_Median_N)[0] + 1
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
                results = []
                for v in range(len(Peak_Begin)):
                    #temp_Final_Peak = pd.DataFrame(columns=['AverageMZ', 'RT', 'RI', 'Int','RTList','SimilarityScore'])
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
                        Int_List_Gaussians = [np.array([gaussian(x, b=RTList[RTPlace[0]], c=c) for x in RTList])for c in gaussian_params]
                        SimilarityScores = [EazyMZDataProcess.CosineSimilarity(Peak_Int_List, g)for g in Int_List_Gaussians]
                        SimilarityScore = max(SimilarityScores)
                        Same_list = list(filter(lambda x:abs(Final_Peak_Detect['AverageMZ'][x]-MZ)/MZ<self.__param['MS1_Tor'] and abs(Final_Peak_Detect['RT'][x]-RT)<2,range(len(Final_Peak_Detect))))
                        if len(Same_list) == 0:
                            temp_Final_Peak = pd.DataFrame([(MZ, RT, RI, Int,RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'RI', 'Int','RTList','SimilarityScore'])
                            results.append(temp_Final_Peak)
                if len(results)>0:
                    temp_Final_Peak = pd.concat(results)
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
                        temp_Final_Peak.sort_values('Int',inplace=True,ignore_index=True)
                        temp_Final_Peak = temp_Final_Peak.iloc[range(Heuristic_peak_list.at[ii,'PeakNumber']),:].copy()
                    Final_Peak_Detect = pd.concat([Final_Peak_Detect, temp_Final_Peak],ignore_index=True)
        return [Final_Peak_Detect.copy(),Heuristic_peak_list.copy()]

    def PeakDetecter(self,UnprocessList,ListIndex,WhetherAppend=True):
        # 单个色谱峰识别模块，由下方另一函数调用
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        Auto_RT_List = UnprocessList['RTList'][ListIndex]
        Auto_Int_List = UnprocessList['Int'][ListIndex]
        if Auto_RT_List[-1]-Auto_RT_List[0] > self.get_param('min_RT_width')*0.8:
            if np.max(Auto_Int_List) > self.get_param('min_Int'):
                Diff_List = np.diff(Auto_Int_List)
                if Diff_List.max() > self.get_param('min_Int')*0.1 and Diff_List.min() < self.get_param('min_Int')*(-0.1):
                    if WhetherAppend==True:
                        RTL = Auto_RT_List[0]-self.get_param('RT_Tor')
                        RTR = Auto_RT_List[-1]+self.get_param('RT_Tor')
                    else:
                        RTL = Auto_RT_List[0]
                        RTR = Auto_RT_List[-1]
                    MZ = np.around(UnprocessList['AverageMZ'][ListIndex],5)
                    [Auto_RT_List, Origin_Int_List] = EazyMZDataProcess.ExtractDataPoint(self,MZ, RTR, RTL,s_min=False, plt_for_test=False, smooth_index=0)
                    if self.get_param('smooth')>0:
                        smooth_index =(self.get_param('smooth')-1)//2
                        Auto_Int_List = list(map(lambda x:GaussSmooth(Origin_Int_List[x-smooth_index:x+smooth_index+1]) if smooth_index<=x<=len(Origin_Int_List)-smooth_index else Origin_Int_List[x],range(len(Origin_Int_List))))
                    else:
                        Auto_Int_List = Origin_Int_List
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
                    #SD_Place = list(filter(lambda x: SD_List[x] < SD_Median, range(len(SD_List))))
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
        
    def detect_Peak(self):
        def gaussian(x,a=1,b=10,c=2):
            y=a*math.exp((-1*(x-b)**2)/(2*c**2))
            return y
        def vectorized_ROI(temp_low, temp_high, keys, values):
            mask = (keys >= temp_low) & (keys <= temp_high)
            selected_values = values[mask]
            return [item for sublist in selected_values for item in sublist]
        # 色谱峰识别
        Final_Peak_Detect = pd.DataFrame(columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
        #Auto_List = pd.DataFrame(columns=['AverageMZ', 'MZ', 'Int', 'RTList'])
        Auto_Peak_Detect = pd.DataFrame(columns=['AverageMZ', 'MZ', 'Int', 'RTList'])
        Origin_MZ_List = self.Origin_MZ_List
        Origin_Int_List = self.Origin_Int_List
        Origin_RT_List = self.Origin_RT_List
        gaussian_params = [0.5, 1, 2, 3 ,4, 5]
        #self.result_PPL = []
        '''------'''
        results = []
        bar = Bar('Processing', max=len(Origin_RT_List))
        for i in range(len(Origin_RT_List)):  
            Potential_Peak_List = pd.DataFrame(columns=['AverageMZ', 'MZ', 'Int', 'RTList'])
            bar.next()
            temp_MZ_List = np.array([])
            temp_Int_List = np.array([])
            for ii in range(len(Origin_MZ_List[i])):
                if Origin_Int_List[i][ii] > self.get_param('max_Noise'):
                    if len(temp_MZ_List) == 0:
                        temp_MZ_List = np.append(temp_MZ_List,Origin_MZ_List[i][ii])
                        temp_Int_List = np.append(temp_Int_List,Origin_Int_List[i][ii])
                    else:
                        if abs(temp_MZ_List[-1]-Origin_MZ_List[i][ii])/temp_MZ_List[-1] > self.get_param('MS1_Tor'):
                            temp_MZ_List = np.append(temp_MZ_List,Origin_MZ_List[i][ii])
                            temp_Int_List = np.append(temp_Int_List,Origin_Int_List[i][ii])
                        else:
                            if temp_Int_List[-1] < Origin_Int_List[i][ii]:
                                temp_MZ_List[-1] = Origin_MZ_List[i][ii]
                                temp_Int_List[-1] = Origin_Int_List[i][ii]
            temp_MZ_List = np.around(temp_MZ_List,5)
            temp_Int_List = np.around(temp_Int_List,1)
            pd.set_option('mode.chained_assignment', None)
            Auto_MZ_List = np.array(Auto_Peak_Detect['AverageMZ'])
            mz_dict = {}
            tolerance = self.get_param('MS1_Tor')
            for idx, mz in enumerate(Auto_MZ_List):
                key = int(mz / tolerance)
                mz_dict.setdefault(key, []).append(idx)
            keys = np.array(list(mz_dict.keys()))
            values = np.array([mz_dict[key] for key in keys])
            temp_ranges = [(iii * (1 - tolerance), iii * (1 + tolerance)) for iii in temp_MZ_List]
            append_Auto_AMZ = []
            append_Auto_MZ = []
            append_Auto_Int = []
            append_Auto_RTList = []
            calculate_mean = []
            df_ROI = pd.DataFrame({'temp_MZ_List':temp_MZ_List,'temp_ranges':temp_ranges,'temp_Int_List':temp_Int_List})
            df_ROI['temp_low'] = np.floor(df_ROI['temp_ranges'].apply(lambda x: x[0] / tolerance).values).astype(int)
            df_ROI['temp_high'] = np.ceil(df_ROI['temp_ranges'].apply(lambda x: x[1] / tolerance).values).astype(int)
            df_ROI['same_MZ_Location'] = [vectorized_ROI(temp_low, temp_high, keys, values) for temp_low, temp_high in zip(df_ROI['temp_low'], df_ROI['temp_high'])]
            temp_MZ_List = df_ROI['temp_MZ_List'].values
            temp_Int_List = df_ROI['temp_Int_List'].values
            same_MZ_Locations = df_ROI['same_MZ_Location'].values
            for iii, (same_MZ_Location, iii_Int, temp_MZ) in enumerate(zip(same_MZ_Locations, temp_Int_List, temp_MZ_List)):
                if len(same_MZ_Location) >= 1:
                    calculate_mean += same_MZ_Location
                    for i_same in same_MZ_Location:
                        # 批量更新 Auto_Peak_Detect
                        if Auto_Peak_Detect.at[i_same, 'RTList'][-1] != Origin_RT_List[i]:
                            Auto_Peak_Detect.at[i_same, 'MZ'].append(temp_MZ)
                            Auto_Peak_Detect.at[i_same, 'Int'].append(iii_Int)
                            Auto_Peak_Detect.at[i_same, 'RTList'].append(Origin_RT_List[i])
                        else:
                            if Auto_Peak_Detect.at[i_same, 'Int'][-1] < iii_Int:
                                Auto_Peak_Detect.at[i_same, 'MZ'][-1] = temp_MZ
                                Auto_Peak_Detect.at[i_same, 'Int'][-1] = iii_Int
                else:
                    # 如果没有匹配，直接追加到新列表
                    append_Auto_AMZ.append(temp_MZ)
                    append_Auto_MZ.append([temp_MZ])
                    append_Auto_Int.append([iii_Int])
                    append_Auto_RTList.append([Origin_RT_List[i]])
            Auto_Peak_Detect.loc[list(set(calculate_mean)),'AverageMZ'] = Auto_Peak_Detect.loc[list(set(calculate_mean)),'MZ'].apply(lambda x:np.mean(x))
            temp_Auto_Detect = pd.DataFrame({'AverageMZ': append_Auto_AMZ, 'MZ': append_Auto_MZ,'Int': append_Auto_Int, 'RTList': append_Auto_RTList})
            temp_Auto_Detect = pd.DataFrame({
            'AverageMZ': append_Auto_AMZ,
            'MZ': append_Auto_MZ,
            'Int': append_Auto_Int,
            'RTList': append_Auto_RTList
            })
            Auto_Peak_Detect = pd.concat([Auto_Peak_Detect, temp_Auto_Detect], ignore_index=True)
            Auto_Peak_Detect.sort_values(by='AverageMZ', ascending=True, inplace=True, ignore_index=True)
            len_mz = Auto_Peak_Detect['MZ'].apply(len)
            mz_values = np.array(Auto_Peak_Detect['MZ'].tolist(),dtype=object)
            avg_mz_values = np.array(Auto_Peak_Detect['AverageMZ'])
            int_values = Auto_Peak_Detect['Int'].apply(np.array).tolist()
            to_drop = []
            for x in range(len(Auto_Peak_Detect) - 1):
                if (
                    len_mz[x] >= 3 and
                    EazyMZDataProcess.ppm_compare(avg_mz_values[x], avg_mz_values[x + 1]) < self.get_param('MS1_Tor') and
                    np.all(mz_values[x][-3:] == mz_values[x + 1][-3:])
                ):
                    if len(int_values[x]) > len(int_values[x + 1]):
                        to_drop.append(x + 1)
                    else:
                        to_drop.append(x)
            Auto_Peak_Detect.drop(index=to_drop, inplace=True)
            Auto_Peak_Detect.reset_index(drop=True, inplace=True)
            ran_Auto_Peak = np.where([rt[-1] != Origin_RT_List[i] for rt in Auto_Peak_Detect['RTList']])[0]
            Auto_Peak_Detect.loc[ran_Auto_Peak, 'Int'] = Auto_Peak_Detect.loc[ran_Auto_Peak, 'Int'].apply(EazyMZDataProcess.add_0)
            Auto_Peak_Detect.loc[ran_Auto_Peak, 'RTList'] = Auto_Peak_Detect.loc[ran_Auto_Peak, 'RTList'].apply(lambda x: EazyMZDataProcess.add_RT(x, Origin_RT_List[i]))
            if len(ran_Auto_Peak)>0:
                temp_PPL = Auto_Peak_Detect.iloc[ran_Auto_Peak].copy()  
                temp_PPL = temp_PPL[temp_PPL['Int'].apply(len) >= self.get_param('Points')*0.5].reset_index(drop=True)
                #self.result_PPL.append(temp_PPL)
                Potential_Peak_List = temp_PPL.copy()
                Potential_Peak_List.reset_index(drop=True,inplace=True)
                for iv in range(len(Potential_Peak_List)):  
                    RT_match = []
                    if len(Potential_Peak_List['Int'][iv])>=self.get_param('Points')/2 and Potential_Peak_List.at[iv,'AverageMZ'] >= self.get_param('min_MZ'):
                        Peak_Result = EazyMZDataProcess.PeakDetecter(self,Potential_Peak_List,iv,WhetherAppend=True)
                        if len(Peak_Result) >0:
                            Peak_Begin = Peak_Result[0]
                            Peak_End = Peak_Result[1]
                            Peak_Top = Peak_Result[2]
                            Auto_Int_List = Peak_Result[3]
                            Auto_RT_List = Peak_Result[4]
                            MZ = Peak_Result[5]
                            for v in range(len(Peak_Begin)):
                                RT = Auto_RT_List[Peak_Top[v]]
                                IntList = Auto_Int_List[Peak_Begin[v]:Peak_End[v]+1]
                                RTList = Auto_RT_List[Peak_Begin[v]:Peak_End[v]+1]
                                Int = np.max(IntList)
                                if RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points') and (IntList[0] < Int*0.20 or IntList[-1] < Int*0.20):
                                    RT_match.append(RT)
                                    Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                    SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                    SimilarityScore = max(SimilarityScores)
                                    temp_Final_Peak = pd.DataFrame([(MZ, RT, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                    results.append(temp_Final_Peak)
                        if len(Peak_Result) > 0 and np.log10(max(Potential_Peak_List['Int'][iv])/self.get_param('min_Int'))>=2:
                            Peak_End.insert(0,0)
                            Peak_Begin.append(len(Auto_RT_List))
                            for v in range(len(Peak_Begin)):
                                if Peak_Begin[v]-Peak_End[v] >= self.get_param('Points') and max(Auto_Int_List[Peak_End[v]:Peak_Begin[v]]) >= self.get_param('min_Int'):
                                    Iterate_Potential_Peak_List = pd.DataFrame({'AverageMZ':MZ,'RTList':[Auto_RT_List[Peak_End[v]:Peak_Begin[v]]],'Int':[Auto_Int_List[Peak_End[v]:Peak_Begin[v]]]})
                                    Peak_Result = EazyMZDataProcess.PeakDetecter(self,Iterate_Potential_Peak_List,0,WhetherAppend=False)
                                    if len(Peak_Result) >0:
                                        Peak_Begin_Iterate = Peak_Result[0]
                                        Peak_End_Iterate = Peak_Result[1]
                                        Peak_Top_Iterate = Peak_Result[2]
                                        Auto_Int_List_Iterate = Peak_Result[3]
                                        Auto_RT_List_Iterate = Peak_Result[4]
                                        MZ_Iterate = Peak_Result[5]
                                        for vi in range(len(Peak_Begin_Iterate)):
                                            RT = Auto_RT_List_Iterate[Peak_Top_Iterate[vi]]
                                            IntList = Auto_Int_List_Iterate[Peak_Begin_Iterate[vi]:Peak_End_Iterate[vi]+1]
                                            RTList = Auto_RT_List_Iterate[Peak_Begin_Iterate[vi]:Peak_End_Iterate[vi]+1]
                                            Int = max(IntList)
                                            if (IntList[0]/Int <= 0.20 or IntList[-1]/Int <= 0.20) and RTList[-1]-RTList[0]>self.get_param('min_RT_width') and len(list(filter(lambda x:abs(RT-x)<self.get_param('RT_Tor')/3,RT_match)))<0:
                                                if RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points'): 
                                                    RT_match.append(RT)
                                                    Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                                    SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                                    SimilarityScore = max(SimilarityScores)
                                                    temp_Final_Peak = pd.DataFrame([(MZ_Iterate, RT, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                                    results.append(temp_Final_Peak)                             
            Auto_Peak_Detect.drop(ran_Auto_Peak,inplace=True)
            Auto_Peak_Detect.reset_index(drop=True,inplace=True)
        bar.finish()
        Potential_Peak_List = Auto_Peak_Detect.copy()
        Potential_Peak_List = Potential_Peak_List[Potential_Peak_List['Int'].apply(len) >= self.get_param('Points')].reset_index(drop=True)
        #self.result_PPL.append(Potential_Peak_List)
        bar = Bar('Peak picking', max=len(Potential_Peak_List))
        for iv in range(len(Potential_Peak_List)):  # len(Auto_Peak_Detect)-1,-1,-1
            RT_match = []
            if len(Potential_Peak_List.at[iv,'Int'])>=self.get_param('Points')/2 and Potential_Peak_List.at[iv,'AverageMZ'] >= self.get_param('min_MZ'):
                Peak_Result = EazyMZDataProcess.PeakDetecter(self,Potential_Peak_List,iv,WhetherAppend=True)
                if len(Peak_Result) >0:
                    Peak_Begin = Peak_Result[0]
                    Peak_End = Peak_Result[1]
                    Peak_Top = Peak_Result[2]
                    Auto_Int_List = Peak_Result[3]
                    Auto_RT_List = Peak_Result[4]
                    MZ = Peak_Result[5]
                    for v in range(len(Peak_Begin)):
                        RT = Auto_RT_List[Peak_Top[v]]
                        IntList = Auto_Int_List[Peak_Begin[v]:Peak_End[v]+1]
                        RTList = Auto_RT_List[Peak_Begin[v]:Peak_End[v]+1]
                        Int = np.max(IntList)
                        if RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points') and (IntList[0] < Auto_Int_List[Peak_Top[v]]*0.20 or IntList[-1] < Auto_Int_List[Peak_Top[v]]*0.20):
                            RT_match.append(RT)
                            Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                            SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                            SimilarityScore = max(SimilarityScores)
                            temp_Final_Peak = pd.DataFrame([(MZ, RT, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                            results.append(temp_Final_Peak)
                if len(Peak_Result) > 0 and np.log10(max(Potential_Peak_List['Int'][iv])/self.get_param('min_Int'))>=2:
                    Peak_End.insert(0,0)
                    Peak_Begin.append(len(Auto_RT_List))
                    for v in range(len(Peak_Begin)):
                        if Peak_Begin[v]-Peak_End[v] >= self.get_param('Points') and max(Auto_Int_List[Peak_End[v]:Peak_Begin[v]]) >= self.get_param('min_Int'):
                            Iterate_Potential_Peak_List = pd.DataFrame({'AverageMZ':MZ,'RTList':[Auto_RT_List[Peak_End[v]:Peak_Begin[v]]],'Int':[Auto_Int_List[Peak_End[v]:Peak_Begin[v]]]})
                            Peak_Result = EazyMZDataProcess.PeakDetecter(self,Iterate_Potential_Peak_List,0,WhetherAppend=False)
                            if len(Peak_Result) >0:
                                Peak_Begin_Iterate = Peak_Result[0]
                                Peak_End_Iterate = Peak_Result[1]
                                Peak_Top_Iterate = Peak_Result[2]
                                Auto_Int_List_Iterate = Peak_Result[3]
                                Auto_RT_List_Iterate = Peak_Result[4]
                                MZ_Iterate = Peak_Result[5]
                                for vi in range(len(Peak_Begin_Iterate)):
                                    RT = Auto_RT_List_Iterate[Peak_Top_Iterate[vi]]
                                    IntList = Auto_Int_List_Iterate[Peak_Begin_Iterate[vi]:Peak_End_Iterate[vi]+1]
                                    RTList = Auto_RT_List_Iterate[Peak_Begin_Iterate[vi]:Peak_End_Iterate[vi]+1]
                                    Int = max(IntList)
                                    if (IntList[0]/Int <= 0.20 or IntList[-1]/Int <= 0.20) and RTList[-1]-RTList[0]>self.get_param('min_RT_width') and len(list(filter(lambda x:abs(RT-x)<self.get_param('RT_Tor')/3,RT_match)))<0:
                                        if RT >= self.get_param('min_RT') and Int > self.get_param('min_Int') and len(RTList)>=self.get_param('Points'): 
                                            RT_match.append(RT)
                                            Int_List_Gaussians = [np.array([gaussian(x, b=RT, c=c) for x in RTList])for c in gaussian_params]
                                            SimilarityScores = [EazyMZDataProcess.CosineSimilarity(IntList, g)for g in Int_List_Gaussians]
                                            SimilarityScore = max(SimilarityScores)
                                            temp_Final_Peak = pd.DataFrame([(MZ_Iterate, RT, Int,[MZ],[Int], IntList, RTList,SimilarityScore)], columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity'])
                                            results.append(temp_Final_Peak)     
        bar.finish()
        #self.result_PPL = pd.concat(self.result_PPL,ignore_index=True)
        #self.result_PPL['RTL'] = self.result_PPL['RTList'].apply(lambda x:x[0])
        #self.result_PPL['RTR'] = self.result_PPL['RTList'].apply(lambda x:x[-1])
        if results:
            Final_Peak_Detect = pd.concat([Final_Peak_Detect]+results,ignore_index=True) 
            data_points = np.column_stack((Final_Peak_Detect['RT'],Final_Peak_Detect['AverageMZ']))
            radius = 1
            kdtree = KDTree(data_points)
            self.del_list = []
            for ii in range(len(Final_Peak_Detect)):   
                Match_List = kdtree.query_ball_point([Final_Peak_Detect.at[ii,'RT'],Final_Peak_Detect.at[ii,'AverageMZ']], r=radius)
                Match_List = list(filter(lambda x:x > ii and abs(Final_Peak_Detect.at[ii,'AverageMZ']-Final_Peak_Detect.at[x,'AverageMZ'])/Final_Peak_Detect.at[ii,'AverageMZ']< self.get_param('MS1_Tor'),Match_List))
                for iii in Match_List:
                    if Final_Peak_Detect.at[ii,'Int'] >= Final_Peak_Detect.at[iii,'Int']:
                        self.del_list.append(iii)
                    else:
                        self.del_list.append(ii)
            self.del_list = list(set(self.del_list))
            Final_Peak_Detect.drop(self.del_list,inplace=True)
            Final_Peak_Detect.reset_index(drop=True,inplace=True)
        self.Final_Peak_Detect = Final_Peak_Detect.copy()
        self.Final_Peak_Detect.sort_values(by='AverageMZ',axis=0,ascending=True,inplace=True,ignore_index=True)
        self.Potential_Peak_List = Potential_Peak_List.copy()
        '''
        self.Auto_List = Auto_List.copy()
        if self.get_param('saveAutoList')==True:
            self.Auto_List['len'] = self.Auto_List['RTList'].apply(lambda x:len(x))
            self.Auto_List['len_sec'] = self.Auto_List['RTList'].apply(lambda x:x[-1]-x[0])
        '''
        if len(self.MS2_Data)>0 and self.get_param('assign MS2') == True:
            self.assign_MS2()
        if 'RIIS' in dir(self):
            self.Calculate_RI()
        return Final_Peak_Detect.copy()
    
    def Deconvolution(self,Ex_Data=[],WhetherDel=False):
        def gaussian(x,a=1,b=10,c=2):
            y=a*math.exp((-1*(x-b)**2)/(2*c**2))
            return y
        ''' de-redundant '''
        #self.DeRedundant()
        ''' Deconvolution '''
        self.DeconvolutionList = pd.DataFrame(columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList'])
        self.Final_Peak_Detect['Group'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:[])
        self.Final_Peak_Detect['Note'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:'')
        self.Final_Peak_Detect['Note Level'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:'')
        ISO_Dict = {100:0.1514,200:0.2271,300:0.3028,400:0.3786,500:0.4543,600:0.53,700:0.6165,800:0.6922,900:0.7679}
        Add_Dict_P = {'Na':21.9819,'NH4':17.0265,'K':37.9559,'CH3OH+H':32.0262,'ACN+H':41.0246,'ACN+Na':63.0085,'2*ACN+H':82.0531} #'2Na-H':3.9639,'IsoProp+H':60.0575,'2*K-H':75.9118,'IsoProp+Na+H':83.0473
        Add_Dict_N = {'Cl':35.9767,'HCOO-':43.9898,'Br':79.9261,} #,'CH3COONa-H':82.0030,'+Na-2H':21.9819,'+K-2H':37.9559,'Hac-H':60.0211,'TFA-H':113.9928,'HNO3-H':62.9957,'HCOONa-H':67.9874
        Dimer_Dict_P = {'2M+H':1.0078,'2M+Na':22.9898,'2M+NH4':18.0344} #'2M+K':38.9637,'2M+ACN+H':42.0343,'2M+ACN+Na':64.1063
        Dimer_Dict_N = {'2M-H':-1.0078,'2M+FA-H':44.9977} #'2M+Hac-H':59.0133
        self.Final_Peak_Detect.sort_values(by='AverageMZ',ascending=True,inplace=True,ignore_index=True)
        bar = Bar('Isotope and adduct', max=len(self.Final_Peak_Detect))
        data_points = np.array([(value, i) for i, value in enumerate(self.Final_Peak_Detect['RT'])])
        radius = self.get_param('RT_Tor')/3
        kdtree = KDTree(data_points[:, 0].reshape(-1, 1))
        self.Final_Peak_Detect['ID'] = list(range(1,len(self.Final_Peak_Detect)+1))
        '''
        for i in range(len(self.Final_Peak_Detect)):
            self.Final_Peak_Detect.at[i,'ID'] = i+1
        '''
        for i in range(len(self.Final_Peak_Detect)):
            bar.next()
            #if self.Final_Peak_Detect.at[i,'Note'] == '':
            RT_Match_List = kdtree.query_ball_point([self.Final_Peak_Detect.at[i,'RT']], r=radius)
            RT_Match_List = list(filter(lambda x:x>i,RT_Match_List))
            #RT_Match_List = list(filter(lambda x:abs(self.Final_Peak_Detect.at[i,'RT']-self.Final_Peak_Detect.at[x,'RT'])<self.get_param('RT_Tor')/3,range(i+1,len(self.Final_Peak_Detect))))
            MZ_Diff_List = list(map(lambda x:self.Final_Peak_Detect.at[x,'AverageMZ']-self.Final_Peak_Detect.at[i,'AverageMZ'],RT_Match_List))
            if self.get_param('Polarity') == 'Positive':
                Dimer_Diff_List = list(map(lambda x:self.Final_Peak_Detect.at[x,'AverageMZ']-2*(self.Final_Peak_Detect.at[i,'AverageMZ']-1.0078),RT_Match_List))
                Multi_Diff_List = list(map(lambda x:(self.Final_Peak_Detect.at[x,'AverageMZ']+1.0078)/2-self.Final_Peak_Detect.at[i,'AverageMZ'],RT_Match_List))
            elif self.get_param('Polarity') == 'Negative':
                Dimer_Diff_List = list(map(lambda x:self.Final_Peak_Detect.at[x,'AverageMZ']-2*(self.Final_Peak_Detect.at[i,'AverageMZ']+1.0078),RT_Match_List))
                Multi_Diff_List = list(map(lambda x:(self.Final_Peak_Detect.at[x,'AverageMZ']-1.0078)/2-self.Final_Peak_Detect.at[i,'AverageMZ'],RT_Match_List))
            for i_Diff in range(len(RT_Match_List)):
                Sim_RTL = min(self.Final_Peak_Detect.at[i,'RTList'][0],self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'RTList'][0])
                Sim_RTR = max(self.Final_Peak_Detect.at[i,'RTList'][-1],self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'RTList'][-1])
                [Origin_RT_List,Origin_Int_List] = self.ExtractDataPoint(self.Final_Peak_Detect.loc[i,'AverageMZ'],Sim_RTR,Sim_RTL,smooth_index=0)
                [Target_RT_List,Target_Int_List] = self.ExtractDataPoint(self.Final_Peak_Detect.loc[RT_Match_List[i_Diff],'AverageMZ'],Sim_RTR,Sim_RTL,smooth_index=0)
                SimilarityScore = EazyMZDataProcess.AdjustedCosineSimilarity(Origin_Int_List,Target_Int_List)
                if SimilarityScore > self.get_param('DeconvolutionSimilarityScore'):
                    #MZ_Diff_Match_ISO = list(filter(lambda x:abs(MZ_Diff_List[i_Diff]-ISO_Dict[x])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor'),list(ISO_Dict.keys())))     
                    Iso_threash = (self.Final_Peak_Detect.at[i,'AverageMZ']//100)*100
                    if Iso_threash <= 900:
                        Iso_threash = ISO_Dict[Iso_threash]
                    if (MZ_Diff_List[i_Diff]%1.0034)/self.Final_Peak_Detect.at[i,'AverageMZ']< self.get_param('MS1_Tor') and MZ_Diff_List[i_Diff]//1.0034 <= 3 and self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int']/self.Final_Peak_Detect.at[i,'Int'] <= Iso_threash:
                        self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = str(int(MZ_Diff_List[i_Diff]//1.0034))+'*C13 -> ' + str(i+1)
                        self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] = 'Isotope'
                        continue
                    elif abs(MZ_Diff_List[i_Diff]%1.0034-1.0034)/self.Final_Peak_Detect.at[i,'AverageMZ']< self.get_param('MS1_Tor') and MZ_Diff_List[i_Diff]//1.0034+1 <= 3 and self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int'] < self.Final_Peak_Detect.at[i,'Int']:
                        self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = str(int(MZ_Diff_List[i_Diff]//1.0034+1))+'*C13 -> ' + str(i+1)
                        self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] = 'Isotope'
                        continue
                    if self.get_param('Polarity') == 'Positive':
                        MZ_Diff_Match_Add = list(filter(lambda x:abs(MZ_Diff_List[i_Diff]-Add_Dict_P[x])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor'),list(Add_Dict_P.keys())))
                        if len(MZ_Diff_Match_Add)>0 and self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int'] < self.Final_Peak_Detect.at[i,'Int']:
                            if self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] == '' and self.Final_Peak_Detect.at[i,'Note Level'] != 'Isotope':
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] = 'Adduct'
                            '''
                            else:
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note']+ '; '+ MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                            '''
                            continue
                        if abs(Multi_Diff_List[i_Diff])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor') and self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int']*2>self.Final_Peak_Detect.at[i,'Int']:
                            if self.Final_Peak_Detect.at[i,'Note Level'] == '':
                                self.Final_Peak_Detect.at[i,'Note'] = '[M+2H] -> ' + str(self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'ID'])
                                self.Final_Peak_Detect.at[i,'Note Level'] = 'Multiply-charged'
                                continue
                        MZ_Diff_Match_Add = list(filter(lambda x:abs(Dimer_Diff_List[i_Diff]-Dimer_Dict_P[x])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor'),list(Dimer_Dict_P.keys())))
                        if len(MZ_Diff_Match_Add)>0 and self.Final_Peak_Detect.at[i,'Int']/10 > self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int']:
                            if self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] == '':
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] = 'Dimer'
                            else:
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note']+ '; '+ MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                            continue
                    elif self.get_param('Polarity') == 'Negative':
                        MZ_Diff_Match_Add = list(filter(lambda x:abs(MZ_Diff_List[i_Diff]-Add_Dict_N[x])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor'),list(Add_Dict_N.keys())))
                        if len(MZ_Diff_Match_Add)>0 and self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int'] < self.Final_Peak_Detect.at[i,'Int']:
                            if self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] == '' and self.Final_Peak_Detect.at[i,'Note Level'] != 'Isotope' :
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] = 'Adduct'
                            '''
                            else:
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note']+'; '+MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                            '''
                            continue    
                        if abs(Multi_Diff_List[i_Diff])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor') and self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int']*2>self.Final_Peak_Detect.at[i,'Int']:
                            if self.Final_Peak_Detect.at[i,'Note Level'] == '':
                                self.Final_Peak_Detect.at[i,'Note'] = '[M-2H] -> ' + str(self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'ID'])
                                self.Final_Peak_Detect.at[i,'Note Level'] = 'Multiply-charged'
                                continue
                        MZ_Diff_Match_Add = list(filter(lambda x:abs(Dimer_Diff_List[i_Diff]-Dimer_Dict_N[x])/self.Final_Peak_Detect.loc[i,'AverageMZ']<self.get_param('MS1_Tor'),list(Dimer_Dict_N.keys())))
                        if len(MZ_Diff_Match_Add)>0 and self.Final_Peak_Detect.at[i,'Int']/10 > self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Int']:
                            if self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] == '':
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note Level'] = 'Dimer'
                            else:
                                self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note'] = self.Final_Peak_Detect.at[RT_Match_List[i_Diff],'Note']+ '; '+ MZ_Diff_Match_Add[0]+' -> ' + str(i+1)
                            continue
        bar.finish()
        self.Final_Peak_Detect.sort_values(by='AverageMZ',ascending=False,inplace=True,ignore_index=True)
        #self.Final_Peak_Detect.reset_index(drop=True,inplace=True)
        bar = Bar('Cluster peak', max=len(self.Final_Peak_Detect))
        for i in range(len(self.Final_Peak_Detect)):
            bar.next()
            if self.Final_Peak_Detect.at[i,'Note'] == '' and self.Final_Peak_Detect.at[i,'RT'] >= self.get_param('min_RT'):
                if len(self.DeconvolutionList) == 0:
                    temp_DeconvolutionList = pd.DataFrame([(len(self.DeconvolutionList)+1,self.Final_Peak_Detect.at[i,'RTList'],self.Final_Peak_Detect.at[i,'AverageMZ'],self.Final_Peak_Detect.at[i,'RT'],self.Final_Peak_Detect.at[i,'Int'],[self.Final_Peak_Detect.at[i,'AverageMZ']],[1],[i])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList'])
                    self.DeconvolutionList = pd.concat([self.DeconvolutionList,temp_DeconvolutionList],ignore_index=True)
                else:
                    Match_List = list(filter(lambda x:abs(self.Final_Peak_Detect.at[i,'RT']-self.DeconvolutionList.at[x,'RT'])<self.get_param('RT_Tor')/3 and self.Final_Peak_Detect.at[i,'Int']/self.DeconvolutionList.at[x,'Int']<0.9 and self.DeconvolutionList.at[x,'MZ']-self.Final_Peak_Detect.at[i,'AverageMZ']>15.9,range(len(self.DeconvolutionList))))
                    if len(Match_List) == 0:
                        temp_DeconvolutionList = pd.DataFrame([(len(self.DeconvolutionList)+1,self.Final_Peak_Detect.at[i,'RTList'],self.Final_Peak_Detect.at[i,'AverageMZ'],self.Final_Peak_Detect.at[i,'RT'],self.Final_Peak_Detect.at[i,'Int'],[self.Final_Peak_Detect.at[i,'AverageMZ']],[1],[i])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList'])
                        self.DeconvolutionList = pd.concat([self.DeconvolutionList,temp_DeconvolutionList],ignore_index=True)
                    else:
                        if_add = False
                        SS_List = []
                        for ii in Match_List:
                            RTL = min(self.DeconvolutionList.at[ii,'RT_List'][0],self.Final_Peak_Detect.at[i,'RTList'][0])
                            RTR = max(self.DeconvolutionList.at[ii,'RT_List'][-1],self.Final_Peak_Detect.at[i,'RTList'][-1])
                            [RT_List_D,Int_List_D] = self.ExtractDataPoint(self.DeconvolutionList.at[ii,'MZ'], RTR, RTL,smooth_index=0) 
                            [RT_List_F,Int_List_F] = self.ExtractDataPoint(self.Final_Peak_Detect.at[i,'AverageMZ'], RTR, RTL,smooth_index=0) 
                            SimilarityScore = EazyMZDataProcess.CosineSimilarity(Int_List_D,Int_List_F)
                            if max(Int_List_F)/max(Int_List_D) < 0.9:
                                SS_List.append(SimilarityScore)
                            else:
                                SS_List.append(0)
                        if max(SS_List) > self.get_param('DeconvolutionSimilarityScore'):
                            if_add = True
                            #ii = SS_List.index(max(SS_List))
                            ii = list(filter(lambda x:SS_List[x]>self.get_param('DeconvolutionSimilarityScore'),range(len(SS_List))))
                            ii = list(filter(lambda x:max(self.DeconvolutionList.loc[np.array(Match_List)[ii],'Int'])==self.DeconvolutionList.loc[Match_List[x],'Int'],ii))[0]
                            '''
                            if self.Final_Peak_Detect.at[i,'Int'] > self.DeconvolutionList.at[ii,'Int']:
                                self.DeconvolutionList.at[Match_List[ii],'MZ'] = self.Final_Peak_Detect.at[i,'AverageMZ']
                                self.DeconvolutionList.at[Match_List[ii],'RT'] = self.Final_Peak_Detect.at[i,'RT']
                                self.DeconvolutionList.at[Match_List[ii],'RT_List'] = self.Final_Peak_Detect.at[i,'RTList']
                                self.DeconvolutionList.at[Match_List[ii],'Int'] = self.Final_Peak_Detect.at[i,'Int']
                            '''
                            self.DeconvolutionList.at[Match_List[ii],'MZ_List'].append(self.Final_Peak_Detect.at[i,'AverageMZ'])
                            self.DeconvolutionList.at[Match_List[ii],'Score'].append(max(SS_List))
                            self.DeconvolutionList.at[Match_List[ii],'IndexList'].append(i)
                        if if_add == False:
                            temp_DeconvolutionList = pd.DataFrame([(len(self.DeconvolutionList)+1,self.Final_Peak_Detect.at[i,'RTList'],self.Final_Peak_Detect.at[i,'AverageMZ'],self.Final_Peak_Detect.at[i,'RT'],self.Final_Peak_Detect.at[i,'Int'],[self.Final_Peak_Detect.at[i,'AverageMZ']],[1],[i])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList'])
                            self.DeconvolutionList = pd.concat([self.DeconvolutionList,temp_DeconvolutionList])
                            self.DeconvolutionList.reset_index(drop=True,inplace=True)
                            #self.Final_Peak_Detect.at[i,'Group'] = [len(self.DeconvolutionList)+1]
        bar.finish()
        del_List = list(filter(lambda x:len(self.DeconvolutionList.at[x,'MZ_List'])==1,range(len(self.DeconvolutionList))))
        self.DeconvolutionList.drop(del_List,inplace=True)
        self.DeconvolutionList.reset_index(drop=True,inplace=True)   
        #self.Final_Peak_Detect['Note Level'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:'')
        self.Final_Peak_Detect['Group'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:0)
        self.DeconvolutionList['IndexListMS2'] = self.DeconvolutionList['IndexList'].apply(lambda x:[])
        ''' MS2 '''
        if 'MS2_MZ' not in self.Final_Peak_Detect.keys():
            self.Final_Peak_Detect['MS2_MZ'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:[])
            self.Final_Peak_Detect['MS2_Int'] = self.Final_Peak_Detect['AverageMZ'].apply(lambda x:[])
        for i in range(len(self.DeconvolutionList)):
            IndexListFix = []
            for ii in self.DeconvolutionList.at[i,'IndexList']:
                add_ILF = True
                for iii in self.DeconvolutionList.at[i,'IndexList']:
                    if len(self.Final_Peak_Detect.loc[iii,'MS2_MZ'])>0 and self.Final_Peak_Detect.loc[ii,'AverageMZ'] < self.Final_Peak_Detect.loc[iii,'AverageMZ']:
                        ISF_List = list(filter(lambda x:abs(self.Final_Peak_Detect.loc[ii,'AverageMZ']-self.Final_Peak_Detect.loc[iii,'MS2_MZ'][x])/self.Final_Peak_Detect.loc[ii,'AverageMZ']<self.get_param('MS1_Tor'),range(len(self.Final_Peak_Detect.loc[iii,'MS2_MZ']))))
                        if len(ISF_List) > 0:
                            add_ILF = False
                            if len(self.Final_Peak_Detect.at[ii,'Note']) > 0:
                                self.Final_Peak_Detect.at[ii,'Note'] =self.Final_Peak_Detect.loc[ii,'Note']+ '; ISF -> Peak '+str(self.Final_Peak_Detect.loc[iii,'ID'])
                            else:
                                self.Final_Peak_Detect.at[ii,'Note'] = 'ISF -> Peak '+str(self.Final_Peak_Detect.loc[iii,'ID'])
                                self.Final_Peak_Detect.at[ii,'Note Level'] = 'MS2'
                if add_ILF == True:
                    IndexListFix.append(ii)
                else:
                    self.DeconvolutionList.at[i,'IndexListMS2'].append(ii)
            self.DeconvolutionList.at[i,'IndexList'] = IndexListFix
        ''' Ex Gradient '''
        if Ex_Data != []:
            Ex_Data.set_param('min_Int',self.get_param('min_Int')*0.5)
            Ex_Data.set_param('RT_Tor',self.get_param('RT_Tor'))
            Ex_Data.set_param('Points',self.get_param('Points')*0.8)
            Ex_Data.set_param('RI_Tor',self.get_param('RI_Tor'))
            Ex_Data.set_param('max_Noise',1000)
            Ex_DeconvolutionList = pd.DataFrame(columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList','Ex_Group','Ex_IndexList'])
            bar = Bar('Ex-gradient peak picking', max=len(self.DeconvolutionList))
            for i in range(len(self.DeconvolutionList)):
                bar.next()
                RTL = 0
                RTR = 0
                for ii in self.DeconvolutionList.at[i,'IndexList']:
                    RI = self.Final_Peak_Detect.at[ii,'RI']
                    if RTL == 0:
                        RTL = EazyMZDataProcess.RI_to_RT(RI*(1-self.get_param('RI_Tor')),Ex_Data.RIIS)
                    else:
                        RTL = min(RTL,EazyMZDataProcess.RI_to_RT(RI*(1-self.get_param('RI_Tor')),Ex_Data.RIIS))
                    RTR = max(RTR,EazyMZDataProcess.RI_to_RT(RI*(1+self.get_param('RI_Tor')),Ex_Data.RIIS))
                Ex_PPL = pd.DataFrame(columns=['AverageMZ', 'Int', 'RTList','IndexNumber','Origin_Int'])
                Ex_FPD = pd.DataFrame(columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity','IndexNumber'])
                for ii in self.DeconvolutionList.at[i,'IndexList']:
                    MZ = self.Final_Peak_Detect.at[ii,'AverageMZ']
                    [RT_List,Int_List] = Ex_Data.ExtractDataPoint(MZ, RTR+30, RTL-30,smooth_index=0) 
                    temp_PPL = pd.DataFrame([(MZ,Int_List,RT_List,ii,self.Final_Peak_Detect.at[ii,'Int'])],columns=['AverageMZ', 'Int', 'RTList','IndexNumber','Origin_Int'])
                    if len(temp_PPL)>0:
                        #Ex_PPL = Ex_PPL.dropna(how='all', axis=1)
                        temp_PPL = temp_PPL.dropna(how='all', axis=1)
                        Ex_PPL = pd.concat([Ex_PPL,temp_PPL])
                Ex_PPL.reset_index(drop=True,inplace=True)
                Not_Find_Ex = []
                for iv in range(len(Ex_PPL)):  # len(Auto_Peak_Detect)-1,-1,-1
                    Find = False
                    if 'Peak_Begin' in dir():
                        del Peak_Begin
                    if len(Ex_PPL['Int'][iv])>self.get_param('Points'):
                        Peak_Result = EazyMZDataProcess.PeakDetecter(Ex_Data,Ex_PPL,iv)
                        if len(Peak_Result) >0:
                            Peak_Begin = Peak_Result[0]
                            Peak_End = Peak_Result[1]
                            Peak_Top = Peak_Result[2]
                            Auto_Int_List = Peak_Result[3]
                            Auto_RT_List = Peak_Result[4]
                            MZ = Peak_Result[5]
                            min_index = list(filter(lambda x:Ex_PPL.at[iv,'Origin_Int']*5>=max(Auto_Int_List[Peak_Begin[x]:Peak_End[x]+1])>=Ex_PPL.at[iv,'Origin_Int']*0.1 and max(Auto_Int_List[Peak_Begin[x]:Peak_End[x]+1]) >= Ex_Data.get_param('min_Int'),range(len(Peak_Top))))
                            min_index = list(filter(lambda x:abs(Auto_RT_List[Peak_Top[x]]-0.5*(RTR+RTL))==min(abs(Auto_RT_List[np.array(Peak_Top)[min_index]]-0.5*(RTR+RTL))),min_index))
                            for v in min_index:
                                RT = Auto_RT_List[Peak_Top[v]]
                                IntList = Auto_Int_List[Peak_Begin[v]:Peak_End[v]+1]
                                RTList = Auto_RT_List[Peak_Begin[v]:Peak_End[v]+1]
                                if max(IntList) > self.get_param('min_Int')*0.7 and RTList[-1]-RTList[0]>self.get_param('min_RT')*0.7:
                                    Int_L = IntList[0]
                                    Int_R = IntList[-1]     
                                    if Int_L < Auto_Int_List[Peak_Top[v]]*0.5 or Int_R < Auto_Int_List[Peak_Top[v]]*0.5:
                                        Int = np.max(IntList)
                                        Same_MZ_list = list(filter(lambda x:abs(Ex_FPD['AverageMZ'][x]-MZ)/MZ<0.000010,range(len(Ex_FPD))))
                                        Same_RT_list = list(filter(lambda x:abs(RT-Ex_FPD['RT'][x])<2/2,range(len(Ex_FPD))))
                                        Same_list = list(filter(lambda x:x in Same_MZ_list,Same_RT_list))
                                        if len(Same_list) > 0:
                                                pass
                                        else:  
                                            Find = True
                                            Int_List_Gaussian_1 = list(map(lambda x:gaussian(x,b=RT,c=1),RTList))
                                            Int_List_Gaussian_2 = list(map(lambda x:gaussian(x,b=RT,c=1.3),RTList))
                                            Int_List_Gaussian_3 = list(map(lambda x:gaussian(x,b=RT,c=1.6),RTList))
                                            Int_List_Gaussian_4 = list(map(lambda x:gaussian(x,b=RT,c=2),RTList))
                                            SimilarityScore = max(EazyMZDataProcess.CosineSimilarity(IntList,Int_List_Gaussian_1),
                                                                 EazyMZDataProcess.CosineSimilarity(IntList,Int_List_Gaussian_2),
                                                                 EazyMZDataProcess.CosineSimilarity(IntList,Int_List_Gaussian_3),
                                                                 EazyMZDataProcess.CosineSimilarity(IntList,Int_List_Gaussian_4))
                                            temp_Final_Peak = pd.DataFrame([(MZ, RT, Int,[MZ],[Int], IntList, RTList,SimilarityScore,Ex_PPL.at[iv,'IndexNumber'])], columns=['AverageMZ', 'RT', 'Int','MS_Range','Int_Range','IntList', 'RTList','GaussianSimilarity','IndexNumber'])
                                            #Ex_FPD = Ex_FPD.dropna(how='all', axis=1)
                                            temp_Final_Peak = temp_Final_Peak.dropna(how='all', axis=1)
                                            Ex_FPD = pd.concat([Ex_FPD, temp_Final_Peak])
                                            Ex_FPD.reset_index(drop=True,inplace=True)
                        if Find == False:
                            Not_Find_Ex.append(Ex_PPL.loc[iv,'IndexNumber'])
                if len(Ex_FPD) > 1:
                    for i_Ex in range(len(Ex_FPD)):
                        if i_Ex == 0:
                            temp_DeconvolutionList = []
                            temp_DeconvolutionList = pd.DataFrame([(1,Ex_FPD.at[i_Ex,'RTList'],Ex_FPD.at[i_Ex,'AverageMZ'],Ex_FPD.at[i_Ex,'RT'],Ex_FPD.at[i_Ex,'Int'],[Ex_FPD.at[i_Ex,'AverageMZ']],[1],self.DeconvolutionList.at[i,'IndexList'],len(Ex_DeconvolutionList)+len(temp_DeconvolutionList)+1,[Ex_FPD.at[i_Ex,'IndexNumber']])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList','Ex_Group','Ex_IndexList'])
                        else:
                            Match_List = list(filter(lambda x:abs(Ex_FPD.at[i_Ex,'RT']-temp_DeconvolutionList.at[x,'RT'])<self.get_param('RT_Tor')/3,range(len(temp_DeconvolutionList))))
                            if len(Match_List) == 0:
                                temp_DL = pd.DataFrame([(len(temp_DeconvolutionList)+1,Ex_FPD.at[i_Ex,'RTList'],Ex_FPD.at[i_Ex,'AverageMZ'],Ex_FPD.at[i_Ex,'RT'],Ex_FPD.at[i_Ex,'Int'],[Ex_FPD.at[i_Ex,'AverageMZ']],[1],self.DeconvolutionList.at[i,'IndexList'],len(Ex_DeconvolutionList)+len(temp_DeconvolutionList)+1,[Ex_FPD.at[i_Ex,'IndexNumber']])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList','Ex_Group','Ex_IndexList'])
                                temp_DeconvolutionList = pd.concat([temp_DeconvolutionList,temp_DL])
                                temp_DeconvolutionList.reset_index(drop=True,inplace=True)
                            else:
                                if_add = False
                                #SS_List = []
                                Valified_List = []
                                for ii in Match_List:
                                    RTL = min(temp_DeconvolutionList.at[ii,'RT_List'][0],Ex_FPD.at[i_Ex,'RTList'][0])
                                    RTR = max(temp_DeconvolutionList.at[ii,'RT_List'][-1],Ex_FPD.at[i_Ex,'RTList'][-1])
                                    RT = 0.5*(temp_DeconvolutionList.at[ii,'RT']+Ex_FPD.at[i_Ex,'RT'])
                                    RT_Diff = min(RT-RTL,RTR-RT)
                                    [RT_List_D,Int_List_D] = Ex_Data.ExtractDataPoint(temp_DeconvolutionList.at[ii,'MZ'], RT+RT_Diff, RT-RT_Diff,smooth_index=5) 
                                    [RT_List_F,Int_List_F] = Ex_Data.ExtractDataPoint(Ex_FPD.at[i_Ex,'AverageMZ'], RT+RT_Diff, RT-RT_Diff,smooth_index=5) 
                                    SimilarityScore = EazyMZDataProcess.AdjustedCosineSimilarity(Int_List_D,Int_List_F)
                                    if max(Int_List_F)/max(Int_List_D) <0.9 and SimilarityScore >= self.get_param('DeconvolutionSimilarityScore'):
                                        Valified_List.append(ii)
                                '''
                                        SS_List.append(SimilarityScore)
                                    else:
                                        SS_List.append(0)
                                if max(SS_List) > self.get_param('DeconvolutionSimilarityScore'):
                                '''
                                if len(Valified_List) > 0:
                                    if_add = True
                                    #ii = SS_List.index(max(SS_List))
                                    ii = list(temp_DeconvolutionList.loc[Valified_List,'Int']).index(max(temp_DeconvolutionList.loc[Valified_List,'Int']))
                                    if Ex_FPD.at[i_Ex,'Int']/temp_DeconvolutionList.at[Match_List[ii],'Int'] < 0.9:
                                        temp_DeconvolutionList.at[Valified_List[ii],'Ex_IndexList'].append(Ex_FPD.at[i_Ex,'IndexNumber'])
                                        temp_DeconvolutionList.at[Valified_List[ii],'MZ_List'].append(Ex_FPD.at[i_Ex,'AverageMZ'])
                                    else:
                                        temp_DL = pd.DataFrame([(len(temp_DeconvolutionList)+1,Ex_FPD.at[i_Ex,'RTList'],Ex_FPD.at[i_Ex,'AverageMZ'],Ex_FPD.at[i_Ex,'RT'],Ex_FPD.at[i_Ex,'Int'],[Ex_FPD.at[i_Ex,'AverageMZ']],[1],self.DeconvolutionList.at[i,'IndexList'],len(temp_DeconvolutionList)+1,[Ex_FPD.at[i_Ex,'IndexNumber']])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList','Ex_Group','Ex_IndexList'])
                                        temp_DeconvolutionList = pd.concat([temp_DeconvolutionList,temp_DL])
                                        temp_DeconvolutionList.reset_index(drop=True,inplace=True)
                                if if_add == False:
                                    temp_DL = pd.DataFrame([(len(temp_DeconvolutionList)+1,Ex_FPD.at[i_Ex,'RTList'],Ex_FPD.at[i_Ex,'AverageMZ'],Ex_FPD.at[i_Ex,'RT'],Ex_FPD.at[i_Ex,'Int'],[Ex_FPD.at[i_Ex,'AverageMZ']],[1],self.DeconvolutionList.at[i,'IndexList'],len(temp_DeconvolutionList)+1,[Ex_FPD.at[i_Ex,'IndexNumber']])],columns=['Group', 'RT_List','MZ','RT','Int','MZ_List','Score','IndexList','Ex_Group','Ex_IndexList'])
                                    temp_DeconvolutionList = pd.concat([temp_DeconvolutionList,temp_DL])
                                    temp_DeconvolutionList.reset_index(drop=True,inplace=True)
                    Group_Mark = []
                    index_Valified = []
                    for ii_Ex in range(len(temp_DeconvolutionList)):
                        if len(temp_DeconvolutionList.at[ii_Ex,'Ex_IndexList'])>1:
                            Group_Number = max(self.Final_Peak_Detect['Group'][~self.Final_Peak_Detect['Group'].apply(lambda x:isinstance(x,str))])+1
                            Group_Mark.append(Group_Number)
                            index_Valified = index_Valified + temp_DeconvolutionList.at[ii_Ex,'Ex_IndexList']
                            Max_MZ = max(self.Final_Peak_Detect.loc[temp_DeconvolutionList.at[ii_Ex,'Ex_IndexList'],'AverageMZ'])
                            for iii_Ex in temp_DeconvolutionList.at[ii_Ex,'Ex_IndexList']:
                                self.Final_Peak_Detect.at[iii_Ex,'Group'] = Group_Number
                                self.Final_Peak_Detect.at[iii_Ex,'Note Level'] = 'Ex-Gradient'
                                if self.Final_Peak_Detect.loc[iii_Ex,'AverageMZ']<Max_MZ:
                                    self.Final_Peak_Detect.at[iii_Ex,'Note'] = 'ISF -> Group '+str(Group_Number)
                        else:
                            Group_Number = max(self.Final_Peak_Detect['Group'][~self.Final_Peak_Detect['Group'].apply(lambda x:isinstance(x,str))])+1
                            Group_Mark.append(Group_Number)
                            index_Valified = index_Valified + temp_DeconvolutionList.at[ii_Ex,'Ex_IndexList']
                    index_Valified = index_Valified + Not_Find_Ex
                    if len(index_Valified)>0:
                        Max_MZ = max(self.Final_Peak_Detect.loc[index_Valified,'AverageMZ'])
                    for ii_Ex in Not_Find_Ex:
                        if self.Final_Peak_Detect.loc[ii_Ex,'AverageMZ'] < Max_MZ:
                            self.Final_Peak_Detect.at[ii_Ex,'Note'] = 'ISF -> Group '+str(Group_Mark[0])
                        self.Final_Peak_Detect.at[ii_Ex,'Group'] = Group_Mark[0]
                        self.Final_Peak_Detect.at[ii_Ex,'Note Level'] = 'MS1'
                elif len(Ex_FPD) <= 1:
                    Group_Number = max(self.Final_Peak_Detect['Group'][~self.Final_Peak_Detect['Group'].apply(lambda x:isinstance(x,str))])+1
                    Max_MZ = max(self.Final_Peak_Detect.loc[self.DeconvolutionList.at[i,'IndexList'],'AverageMZ'])
                    for ii in self.DeconvolutionList.at[i,'IndexList']:
                        self.Final_Peak_Detect.at[ii,'Group'] = Group_Number
                        self.Final_Peak_Detect.at[ii,'Note Level'] = 'MS1'
                        if self.Final_Peak_Detect.loc[ii,'AverageMZ']<Max_MZ:
                            self.Final_Peak_Detect.at[ii,'Note'] = 'ISF -> Group '+str(Group_Number)
            bar.finish()
        else:
            for i in range(len(self.DeconvolutionList)):
                for ii in range(len(self.DeconvolutionList.at[i,'IndexList'])):
                    self.Final_Peak_Detect.at[self.DeconvolutionList.at[i,'IndexList'][ii],'Group'] = i+1
                    if ii == 0 :
                        self.Final_Peak_Detect.at[self.DeconvolutionList.at[i,'IndexList'][ii],'Note'] = 'Pre'
                    else:
                        self.Final_Peak_Detect.at[self.DeconvolutionList.at[i,'IndexList'][ii],'Note Level'] = 'MS1'
        self.Final_Peak_Detect.sort_values(by='ID',ascending=True,inplace=True)
        self.Final_Peak_Detect.reset_index(drop=True,inplace=True)
        for i in range(len(self.Final_Peak_Detect)):
            if self.Final_Peak_Detect.at[i,'Note Level'] == 'MS2':
                if self.Final_Peak_Detect.at[i,'Note'].find(';') == -1:
                    end_group = len(self.Final_Peak_Detect.at[i,'Note'])
                else:
                    end_group = self.Final_Peak_Detect.at[i,'Note'].find(';')
                ii = int(self.Final_Peak_Detect.at[i,'Note'][12:end_group])
                if self.Final_Peak_Detect.at[ii-1,'Group'] == 0:
                    self.Final_Peak_Detect.at[ii-1,'Group'] = max(self.Final_Peak_Detect['Group'][~self.Final_Peak_Detect['Group'].apply(lambda x:isinstance(x,str))])+1
                self.Final_Peak_Detect.at[i,'Group'] = self.Final_Peak_Detect.at[ii-1,'Group']
        for i in range(len(self.Final_Peak_Detect)):
            if self.Final_Peak_Detect.at[i,'Group'] == 0:
                self.Final_Peak_Detect.at[i,'Group'] = ''
        self.Final_Peak_Detect.reset_index(drop=True,inplace=True)  
        if WhetherDel == True:
            self.Final_Peak_Detect = self.Final_Peak_Detect[self.Final_Peak_Detect['Note Level'].isin([''])]
            self.Final_Peak_Detect.reset_index(drop=True,inplace=True)  
        
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
    def Find_FContinuous(FD_List, FD_P, Diff_P, FDMode='P' ,FC_Number=2 ,mergeRule='Intersection'):
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
        #self.Final_Peak_Detect['S/N']=self.Final_Peak_Detect['AverageMZ'].apply(lambda x:0)
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
        if type(RIIS) == str:
            self.RIIS = pd.read_excel(RIIS)
            if 'RT' not in list(self.RIIS.keys()):
                self.RIIS['RT']=0.00
                self.RIIS['Int']=0.00
                self.RIIS['Candidate_RT'] = ''
                self.RIIS['Candidate_Int'] = ''
                RTL = 60
                RTR = self.Origin_RT_List[-1]
                for i in range(len(self.RIIS)):                
                    MZ = self.RIIS['m/z'][i]
                    [RT_List, Int_List]=self.ExtractDataPoint(MZ,RTR,RTL)
                    RT_place_1 = Int_List.index(max(Int_List))
                    Int_List_d = Int_List[0:RT_place_1-int(self.__param['Points'])]+Int_List[RT_place_1+int(self.__param['Points']):-1]
                    RT_place_2 = Int_List_d.index(max(Int_List_d))
                    if RT_place_2 >= RT_place_1-int(self.__param['Points']):
                        RT_place_2 += int(self.__param['Points'])*2
                    if Int_List[RT_place_1]>=self.__param['min_Int']:
                        self.RIIS.at[i,'Candidate_RT'] = [RT_List[RT_place_1],RT_List[RT_place_2]]
                        self.RIIS.at[i,'Candidate_Int'] = [Int_List[RT_place_1],Int_List[RT_place_2]]
                        self.RIIS.at[i,'RT'] = RT_List[RT_place_1].copy()
                        self.RIIS.at[i,'Int'] = Int_List[RT_place_1].copy()
                    if Int_List[RT_place_1]/Int_List[RT_place_2] >=5:
                        RTL = RT_List[RT_place_1]
                del_list = list(filter(lambda x:self.RIIS['RT'][x]==0,range(len(self.RIIS))))
                self.RIIS.drop(del_list,inplace=True)
                self.RIIS.reset_index(drop=True,inplace=True)   
                del_fix = 0
                for i in range(len(self.RIIS.iloc[:,0])): #len(self.RIIS.iloc[:,0])
                    i = i - del_fix
                    if i == 0:
                        if self.RIIS['RT'][i] >= self.RIIS['RT'][i+1] and len(self.RIIS['Candidate_RT'][i]) == 2:
                            self.RIIS.at[i,'RT'] = self.RIIS['Candidate_RT'][i][1]
                            if self.RIIS['RT'][i] >= self.RIIS['RT'][i+1]:
                                self.RIIS.drop(i,inplace=True)
                                self.RIIS.reset_index(drop=True,inplace=True)  
                                del_fix = del_fix + 1
                        elif self.RIIS['RT'][i] >= self.RIIS['RT'][i+1] and len(self.RIIS['Candidate_RT'][i]) == 1:
                            self.RIIS.drop(i,inplace=True)
                            self.RIIS.reset_index(drop=True,inplace=True)  
                            del_fix = del_fix + 1
                    elif i == len(self.RIIS.iloc[:,0])-1:
                        if self.RIIS['RT'][i] <= self.RIIS['RT'][i-1] and len(self.RIIS['Candidate_RT'][i]) == 2:
                            self.RIIS.at[i,'RT'] = self.RIIS['Candidate_RT'][i][1]
                            if self.RIIS['RT'][i] <= self.RIIS['RT'][i-1]:
                                self.RIIS.drop(i,inplace=True)
                                self.RIIS.reset_index(drop=True,inplace=True)  
                                del_fix = del_fix + 1
                        elif self.RIIS['RT'][i] <= self.RIIS['RT'][i-1] and len(self.RIIS['Candidate_RT'][i]) == 1:
                            self.RIIS.drop(i,inplace=True)
                            self.RIIS.reset_index(drop=True,inplace=True)  
                            del_fix = del_fix + 1
                    else:
                        if self.RIIS['RT'][i] > self.RIIS['RT'][i-1] and self.RIIS['RT'][i] < self.RIIS['RT'][i+1]:
                            continue
                        if self.RIIS['RT'][i] < self.RIIS['RT'][i-1]:
                            if len(self.RIIS['Candidate_RT'][i]) == 2:
                                self.RIIS.at[i,'RT'] = self.RIIS['Candidate_RT'][i][1]
                                if self.RIIS['RT'][i] <= self.RIIS['RT'][i-1] or self.RIIS['RT'][i] >=self.RIIS['RT'][i+1]:
                                    self.RIIS.drop(i,inplace=True)
                                    self.RIIS.reset_index(drop=True,inplace=True)  
                                    del_fix = del_fix + 1 
                        elif self.RIIS['RT'][i] > self.RIIS['RT'][i+1]:
                            if i >=2 :
                                i_n = ((self.RIIS['RT'][i]-self.RIIS['RT'][i-1])/(self.RIIS['C'][i]-self.RIIS['C'][i-1]))/((self.RIIS['RT'][i-1]-self.RIIS['RT'][i-2])/(self.RIIS['C'][i-1]-self.RIIS['C'][i-2]))
                                i_n_1 = ((self.RIIS['RT'][i+1]-self.RIIS['RT'][i-1])/(self.RIIS['C'][i+1]-self.RIIS['C'][i-1]))/((self.RIIS['RT'][i-1]-self.RIIS['RT'][i-2])/(self.RIIS['C'][i-1]-self.RIIS['C'][i-2]))
                                if abs(1-i_n) > abs(1-i_n_1):
                                    if len(self.RIIS['Candidate_RT'][i]) == 2:
                                            self.RIIS.at[i,'RT'] = self.RIIS['Candidate_RT'][i][1]
                                            if self.RIIS['RT'][i] <= self.RIIS['RT'][i-1] or self.RIIS['RT'][i] >=self.RIIS['RT'][i+1]:
                                                self.RIIS.drop(i,inplace=True)
                                                self.RIIS.reset_index(drop=True,inplace=True)  
                                                del_fix = del_fix + 1 
                            elif len(self.RIIS['Candidate_RT'][i]) == 2:
                                    self.RIIS.at[i,'RT'] = self.RIIS['Candidate_RT'][i][1]
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
                for i in range(len(self.RIIS)):                
                    MZ = self.RIIS['m/z'][i]
                    [RT_List, Int_List]=self.ExtractDataPoint(MZ,self.RIIS['RT'][i]*60+self.get_param('RT_Tor'),self.RIIS['RT'][i]*60-self.get_param('RT_Tor'))
                    self.RIIS['RT'][i] = [RT_List[x] for x in range(len(Int_List)) if Int_List[x]==max(Int_List)][0]
        else:
            self.RIIS = RIIS.copy()
            for i in range(len(self.RIIS)):                
                MZ = self.RIIS['m/z'][i]
                [RT_List, Int_List]=self.ExtractDataPoint(MZ,self.RIIS['RT'][i]*60+self.get_param('RT_Tor'),self.RIIS['RT'][i]*60-self.get_param('RT_Tor'))
                self.RIIS['RT'][i] = [RT_List[x] for x in range(len(Int_List)) if Int_List[x]==max(Int_List)][0]
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
    
    def load_FPD(self,Path):
        with open(Path,'rb') as f:
            self.Final_Peak_Detect = pickle.load(f)
    
    def save_FPD(self,Path):
        with open(Path,'wb') as f:
            pickle.dump(self.Final_Peak_Detect,f)

class DataAlignment(object):
    def __init__(self):
        self.DataBase = pd.DataFrame(columns=['Data','Data_Name','Tag','Final_Peak_Detect'])
        self.AlignmentParam = {'MZ_Tor':0.000010,'RT_Tor':6,'A':0.5,'RI_Tor':0.02,
                               'Miss_Filter':0.8,'Threshold':5,'RI_Alignment':False}
        self.RefList=pd.DataFrame(columns=['m/z','RT','MS_List','RT_List'])
        self.RefList['m/z'] = self.RefList['m/z'].map(lambda x:'%.4f'%x)
    def add_Data_new(self,Data,Tag='Sample'):
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
    def add_Data(self,Data,Tag='Sample'):
        # 添加数据
        if Tag == 'Sample' or  Tag == 'QC':
            def namestr(obj, namespace):
                return [name for name in namespace if namespace[name] is obj]
            self.DataBase.at[len(self.DataBase),'Data']= Data['Data']
            self.DataBase.at[len(self.DataBase)-1,'Data_Name'] = Data['Data'].DataName
            self.DataBase.at[len(self.DataBase)-1,'Tag'] = Tag
            self.DataBase.sort_values(by='Tag',ascending=False,inplace=True,ignore_index=True)
        else:
            print('Data Tag should be Sample or Blank or QC')
    def add_raw_Data(self,DataPath,Tag='Sample',Name=''):
        if Tag == 'Sample' or Tag == 'Blank':
            self.DataBase.loc[len(self.DataBase),'Tag'] = Tag
            self.DataBase.loc[len(self.DataBase)-1,'Data'] = EazyMZDataProcess(DataPath)
            #self.DataBase.loc[len(self.DataBase)-1,'Gradient'] = 0
            if Name == '':
                endName = DataPath.rfind('.')
                beginName = DataPath.rfind('/')
                self.DataBase.loc[len(self.DataBase)-1,'Data_Name'] = DataPath[beginName+1:endName]
            else:
                self.DataBase.loc[len(self.DataBase)-1,'Data_Name'] = Name
            #self.DataBase.loc[len(self.DataBase)-1,'Data'] = EazyMZDataProcess(self.DataBase.loc[i,'Data'])
    def show_Data(self):
        # 显示数据
        print(self.DataBase.iloc[:,[1,2]])
    def get_param(self,key=''):
        if len(key) == 0:
            print(self.AlignmentParam)
        else:
            return self.AlignmentParam[key]
    def set_param(self,Name,Value):
        # 更改参数
        if Name in self.AlignmentParam:
            self.AlignmentParam[Name]=Value
        else:
            print('no such param')
    def RenewRefList(self):
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
            if self.DataBase.at[i,'Tag'] != 'Blank':
                self.RefList.sort_values('m/z',ascending=(True),inplace=True)
                self.RefList.reset_index(drop=True,inplace=True)
                temp_Data = self.DataBase.at[i,'Data'].Final_Peak_Detect.copy()
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
                init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True,inplace=True)
                for ii in range(len(self.RefList)):
                    init_list = range(bisect.bisect_left(init_RefMZ['m/z'],self.RefList.at[ii,'m/z']*(1-self.AlignmentParam['MZ_Tor'])),bisect.bisect_right(init_RefMZ['m/z'],self.RefList.at[ii,'m/z']*(1+self.AlignmentParam['MZ_Tor'])))
                    if len(init_list) == 1:
                        continue
                    elif len(init_list) >1:
                        for iii in init_list:
                            if ii not in init_RefMZ.at[iii,'RowIndex']:
                                init_RefMZ.at[iii,'MZList'].append(self.RefList.at[ii,'m/z'])
                                init_RefMZ.at[iii,'RowIndex'].append(ii)
                                init_RefMZ.at[iii,'m/z']=np.mean(init_RefMZ.at[iii,'MZList'])
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
                    bar = Bar('Alignment '+str(i+1)+' / '+str(len(self.DataBase)), max=len(init_SampleMZ))
                    init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True,inplace=True)
                    for ii in range(len(init_SampleMZ)):
                        bar.next()
                        Sample_MZ = init_SampleMZ.at[ii,'m/z']
                        match_Index = range(bisect.bisect_left(init_RefMZ['m/z'],Sample_MZ*(1-self.AlignmentParam['MZ_Tor']*2)),bisect.bisect_right(init_RefMZ['m/z'],Sample_MZ*(1+self.AlignmentParam['MZ_Tor']*2)))
                        #match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
                        if len(match_Index)>1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                            for i_mI in range(1,len(match_Index)):
                                Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                            Ref_Index = list(set(Ref_Index))
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
                    bar.finish()
                else:
                    Sample_RI = temp_Data.loc[:,'RI']
                    init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True,inplace=True)
                    bar = Bar('Alignment '+str(i+1)+' / '+str(len(self.DataBase)), max=len(init_SampleMZ))
                    for ii in range(len(init_SampleMZ)):
                        bar.next()
                        Sample_MZ = init_SampleMZ.at[ii,'m/z']
                        #match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
                        match_Index = range(bisect.bisect_left(init_RefMZ['m/z'],Sample_MZ*(1-self.AlignmentParam['MZ_Tor']*2)),bisect.bisect_right(init_RefMZ['m/z'],Sample_MZ*(1+self.AlignmentParam['MZ_Tor']*2)))
                        if len(match_Index)>1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                            for i_mI in range(1,len(match_Index)):
                                Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                            Ref_Index = list(set(Ref_Index))
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
                    bar.finish()
                
                add_RefList = pd.DataFrame({'m/z': add_MZ,'RT': add_RT,'RI':add_RI,self.DataBase.at[i, 'Data_Name']: add_Int,
                                            'MS_List':add_MS_List,'RT_List':add_RT_List,'RI_List':add_RI_List})
                self.RefList = pd.concat([self.RefList, add_RefList])
                self.RefList.reset_index(drop=True, inplace=True)
                self.RefList = self.RefList.fillna(0)
        Sample_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Sample'])
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        self.RefList['max_SampleInt'] = self.RefList.apply(lambda x:max(x[Sample_Name]),axis=1)
        self.RefList['Int'] = self.RefList.apply(lambda x:np.mean(x[Sample_Name]),axis=1)
        self.RefList['max_BlankInt'] = 0
        if len(Blank_Name) > 0:
            self.RefList['max_BlankInt'] = self.RefList.apply(lambda x:max(x[Blank_Name]) if max(x[Blank_Name])>0 else 1,axis=1)
    def RenewRefList_new(self):
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
                init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True,inplace=True)
                for ii in range(len(self.RefList)):
                    init_list = range(bisect.bisect_left(init_RefMZ['m/z'],self.RefList.at[ii,'m/z']*(1-self.AlignmentParam['MZ_Tor'])),bisect.bisect_right(init_RefMZ['m/z'],self.RefList.at[ii,'m/z']*(1+self.AlignmentParam['MZ_Tor'])))
                    if len(init_list) == 1:
                        continue
                    elif len(init_list) >1:
                        for iii in init_list:
                            if ii not in init_RefMZ.at[iii,'RowIndex']:
                                init_RefMZ.at[iii,'MZList'].append(self.RefList.at[ii,'m/z'])
                                init_RefMZ.at[iii,'RowIndex'].append(ii)
                                init_RefMZ.at[iii,'m/z']=np.mean(init_RefMZ.at[iii,'MZList'])
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
                    bar = Bar('Alignment '+str(i+1)+' / '+str(len(self.DataBase)), max=len(init_SampleMZ))
                    init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True,inplace=True)
                    for ii in range(len(init_SampleMZ)):
                        bar.next()
                        Sample_MZ = init_SampleMZ.at[ii,'m/z']
                        #match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
                        match_Index = range(bisect.bisect_left(init_RefMZ['m/z'],Sample_MZ*(1-self.AlignmentParam['MZ_Tor']*2)),bisect.bisect_right(init_RefMZ['m/z'],Sample_MZ*(1+self.AlignmentParam['MZ_Tor']*2)))
                        if len(match_Index)>1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                            for i_mI in range(1,len(match_Index)):
                                Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                            Ref_Index = list(set(Ref_Index))
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
                    bar.finish()
                else:
                    Sample_RI = temp_Data.loc[:,'RI']
                    bar = Bar('Alignment '+str(i+1)+' / '+str(len(self.DataBase)), max=len(init_SampleMZ))
                    init_RefMZ.sort_values('m/z',ascending=True,ignore_index=True,inplace=True)
                    for ii in range(len(init_SampleMZ)):
                        bar.next()
                        Sample_MZ = init_SampleMZ.at[ii,'m/z']
                        #match_Index = list(filter(lambda x:abs(init_RefMZ.at[x,'m/z']-Sample_MZ)/Sample_MZ<self.AlignmentParam['MZ_Tor'],range(len(init_RefMZ))))
                        match_Index = range(bisect.bisect_left(init_RefMZ['m/z'],Sample_MZ*(1-self.AlignmentParam['MZ_Tor'])),bisect.bisect_right(init_RefMZ['m/z'],Sample_MZ*(1+self.AlignmentParam['MZ_Tor'])))
                        if len(match_Index)>1:
                            Ref_Index = init_RefMZ.at[match_Index[0],'RowIndex']
                            for i_mI in range(1,len(match_Index)):
                                Ref_Index = Ref_Index + init_RefMZ.at[match_Index[i_mI],'RowIndex']
                            Ref_Index = list(set(Ref_Index))
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
                    bar.finish()
                
                add_RefList = pd.DataFrame({'m/z': add_MZ,'RT': add_RT,'RI':add_RI,self.DataBase.at[i, 'Data_Name']: add_Int,
                                            'MS_List':add_MS_List,'RT_List':add_RT_List,'RI_List':add_RI_List})
                self.RefList = pd.concat([self.RefList, add_RefList])
                self.RefList.reset_index(drop=True, inplace=True)
                self.RefList = self.RefList.fillna(0)
        Sample_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Sample'])
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        self.RefList['max_SampleInt'] = self.RefList.apply(lambda x:max(x[Sample_Name]),axis=1)
        self.RefList['Int'] = self.RefList.apply(lambda x:np.mean(x[Sample_Name]),axis=1)
        self.RefList['max_BlankInt'] = 0
        if len(Blank_Name) > 0:
            self.RefList['max_BlankInt'] = self.RefList.apply(lambda x:max(x[Blank_Name]) if max(x[Blank_Name])>0 else 1,axis=1)
    def BlankFilter(self):
        # 使用Blank数据进行过滤
        Blank_Name = list(self.DataBase['Data_Name'][self.DataBase['Tag']=='Blank'])
        if len(Blank_Name) > 0:
            self.RefList.drop(self.RefList[self.RefList['mix_SampleInt']/self.RefList['mix_BlankInt']<self.AlignmentParam['Threshold']].index,inplace=True)
            self.RefList.reset_index(drop=True,inplace=True)
    def Filter_MissingValue(self,WhetherDel=True):
        # 缺失值删除
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
    def Output_mgf(self,OutputFilePath):
        mgf_output = ''
        csv_mgf = pd.DataFrame(columns=['row ID','row m/z','row retention time','correlation group ID','annotation network number','best ion','auto MS2 verify','identified by n=','partners','neutral M mass','Peak height'])
        for i in  range(len(self.RefList)):
            MZ = self.RefList['m/z'][i]
            RT = np.around(self.RefList['RT'][i]/60,3)
            Int = self.RefList['Int'][i]
            temp_csv_mgf = pd.DataFrame([(i+1,MZ,RT,'','','','','','','',Int)],columns=['row ID','row m/z','row retention time','correlation group ID','annotation network number','best ion','auto MS2 verify','identified by n=','partners','neutral M mass','Peak height'])
            csv_mgf = pd.concat([csv_mgf,temp_csv_mgf])
            mgf_output = mgf_output+'BEGIN IONS\nFEATURE_ID='+str(i+1)+'\nPEPMASS='+str(MZ)+'\nSCANS='+str(i+1)+'\nRTINSECONDS='+str(RT)+'\nCHARGE=1+\nMSLEVEL=2\n'
            for ii in range(len(self.RefList.at[i,'MS2_MZ'])):
                mgf_output = mgf_output+str(self.RefList.at[i,'MS2_MZ'][ii])+' '+str(self.RefList.at[i,'MS2_Int'][ii])+'\n'
            mgf_output = mgf_output+'END IONS\n\n'
        def namestr(obj, namespace):
            return [name for name in namespace if namespace[name] is obj]    
        with open(OutputFilePath + namestr(self,globals())[0]+'.mgf','w')as mgfFile:
            mgfFile.write(mgf_output)
        csv_mgf.to_csv(OutputFilePath + namestr(self,globals())[0]+'.csv',index=False)
    def Output_Result(self,OutputPath):
        self.RefList.to_excel(OutputPath,index=False)
    def RI_Correct(self,RIIS_Path='',RefNumber=0):
        def RTO_to_RTC(RTO,t_n_S,t_n1_S,t_n_R,t_n1_R):
            RTC = (RTO-t_n_S)/(t_n1_S-t_n_S)*(t_n1_R-t_n_R)+t_n_R
            return RTC
        if hasattr(self.DataBase.at[RefNumber,'Data'],'RIIS') == False:
            self.RIIS = pd.read_excel(RIIS_Path)
        else:
            self.RIIS = self.DataBase.at[RefNumber,'Data'].RIIS.copy()
        bar = Bar('Calculate RI',max = len(self.DataBase))
        for i in range(len(self.DataBase)):
            bar.next()
            if hasattr(self.DataBase.at[i,'Data'],'RIIS') == False:
                self.DataBase.at[i,'Data'].RIIS = self.RIIS.copy()
                self.DataBase.at[i,'Data'].RIIS['RT']=0.00
                RTL = self.DataBase.at[i,'Data'].Origin_RT_List[0]
                RTR = self.DataBase.at[i,'Data'].Origin_RT_List[-1]
                for ii in range(len(self.DataBase.at[i,'Data'].RIIS)):                
                    MZ = self.DataBase.at[i,'Data'].RIIS['m/z'][ii]
                    [RT_List, Int_List]=self.DataBase.at[i,'Data'].ExtractDataPoint(MZ,RTR,RTL)
                    RT_place = np.where(Int_List==max(Int_List))[0]
                    self.DataBase.at[i,'Data'].RIIS['RT'][ii] = RT_List[RT_place[0]]
        bar.finish()
        self.RIIS = self.DataBase.at[RefNumber,'Data'].RIIS.copy()
        bar = Bar('Calculate RI',max = len(self.DataBase))
        for i in range(len(self.DataBase)):
            bar.next()
            if i != RefNumber:
                exp = pyopenms.MSExperiment()
                for ii in range(len(self.DataBase.at[i,'Data'].Origin_RT_List)):   
                    spectrum = self.DataBase.at[i,'Data'].OriginData[ii]
                    RTO = self.DataBase.at[i,'Data'].Origin_RT_List[ii]
                    n_place = np.where(self.DataBase.at[i,'Data'].RIIS['RT']<=RTO)[0]
                    n1_place = np.where(self.DataBase.at[i,'Data'].RIIS['RT']>RTO)[0]
                    if len(n_place)==0:
                        n_place = n1_place[0]
                        n1_place = n1_place[1]
                    elif len(n1_place)==0:
                        n1_place = n_place[-1]
                        n_place = n_place[-2] 
                    else:
                        n_place = n_place[-1]
                        n1_place = n1_place[0]
                    t_n_S = self.DataBase.at[i,'Data'].RIIS['RT'][n_place]
                    t_n1_S = self.DataBase.at[i,'Data'].RIIS['RT'][n1_place]
                    t_n_R = self.RIIS['RT'][n_place]
                    t_n1_R = self.RIIS['RT'][n1_place]
                    RTC = RTO_to_RTC(RTO,t_n_S,t_n1_S,t_n_R,t_n1_R)
                    spectrum.setRT(RTC)
                    exp.addSpectrum(spectrum)
                dot_place = self.DataBase.at[i,'Data'].file_path.rfind('.')
                OutputPath = self.DataBase.at[i,'Data'].file_path[0:dot_place]
                pyopenms.MzMLFile().store(OutputPath+'_RI.mzML',exp)
        bar.finish()
    def MultipleGradientReference(self,RefNumber=0):
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
        self.MGRef = pd.DataFrame(columns=['m/z','Int','RI','MS_List','RT_List','MS2_Int','MS2_MZ'])
        for i in range(len(self.DataBase)):
            if self.DataBase.at[i,'Tag'] =='QC':
                self.MGRef[self.DataBase.at[i,'Data_Name']] = 0
        for i in range(len(self.DataBase)):
            if self.DataBase.at[i,'Tag'] =='QC':
                self.MGRef[self.DataBase.at[i,'Data_Name']+'Score'] = 0
        for i in range(len(self.DataBase)):
            if self.DataBase.at[i,'Tag'] =='QC':
                bar = Bar('Alignment '+str(i+1)+' / '+str(len(self.DataBase)), max=len(self.DataBase.at[i,'Data'].Final_Peak_Detect))
                add_MZ = []
                add_Int = []
                add_RT = []
                add_Score = []
                add_MS2_Int = []
                add_MS2_MZ = []
                add_MS_List = []
                add_RT_List = []
                for ii in range(len(self.DataBase.at[i,'Data'].Final_Peak_Detect)):
                    bar.next()
                    Ref_MZ = self.DataBase.at[i,'Data'].Final_Peak_Detect.at[ii,'AverageMZ']
                    Ref_RT = self.DataBase.at[i,'Data'].Final_Peak_Detect.at[ii,'RI']
                    Ref_Int = self.DataBase.at[i,'Data'].Final_Peak_Detect.at[ii,'Int']
                    add_index = list(filter(lambda x:abs(self.MGRef.at[x,'m/z']-Ref_MZ)/Ref_MZ<self.AlignmentParam['MZ_Tor'] and abs(self.MGRef.at[x,'RI']-Ref_RT)/Ref_RT<self.AlignmentParam['RI_Tor'],range(len(self.MGRef))))
                    #add_index = list(filter(lambda x:abs(self.MGRef.at[x,'m/z']-Ref_MZ)/Ref_MZ<self.AlignmentParam['MZ_Tor'] and abs(self.MGRef.at[x,'RI']-Ref_RT)<40,range(len(self.MGRef))))
                    if len(add_index)==0:
                        add_MZ.append(Ref_MZ)
                        add_Int.append(Ref_Int)
                        add_RT.append(Ref_RT)
                        add_Score.append(1)
                        add_MS_List.append([Ref_MZ])
                        add_RT_List.append([Ref_RT])
                        if self.DataBase.at[i,'Tag']=='Sample' and 'MS2_Int' in self.DataBase.at[i,'Data'].Final_Peak_Detect.keys():
                            if len(self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_Int'])>0:
                                add_MS2_Int.append(self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_Int']/max(self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_Int']))
                                add_MS2_MZ.append(self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_MZ'])
                            else:
                                add_MS2_Int.append([])
                                add_MS2_MZ.append([])
                        else:
                            add_MS2_Int.append([])
                            add_MS2_MZ.append([])
                      
                    elif len(add_index)>0:
                        FS_Score = []
                        for iii in range(len(add_index)):
                            FS_Score.append(self.AlignmentParam['A']*np.exp(-0.5*((Ref_RT-self.MGRef.at[add_index[iii],'RI'])/(Ref_RT*self.AlignmentParam['RI_Tor']))**2)+(1-self.AlignmentParam['A'])*np.exp(-0.5*((Ref_MZ-self.MGRef.at[add_index[iii],'m/z'])/(Ref_MZ*self.AlignmentParam['MZ_Tor']))**2))
                        match_table = pd.DataFrame({'add_index':add_index,'FS_Score':FS_Score})   
                        match_table.sort_values('FS_Score',ascending=(False),inplace=True)
                        match_table.reset_index(inplace=True)
                        #takeplace = True
                        if self.MGRef.at[match_table.at[0, 'add_index'], self.DataBase.at[i, 'Data_Name']+'Score'] < match_table.at[0, 'FS_Score']:
                            self.MGRef.at[match_table.at[0, 'add_index'],self.DataBase.at[i, 'Data_Name']] = Ref_Int
                            self.MGRef.at[match_table.at[0, 'add_index'],self.DataBase.at[i, 'Data_Name']+'Score'] = match_table.at[0, 'FS_Score']
                            self.MGRef.at[match_table.at[0, 'add_index'],'MS_List'].append(Ref_MZ)
                            self.MGRef.at[match_table.at[0, 'add_index'],'RT_List'].append(Ref_RT)
                            self.MGRef.at[match_table.at[0, 'add_index'],'m/z'] = np.mean(self.MGRef.at[match_table.at[0, 'add_index'],'MS_List'])
                            self.MGRef.at[match_table.at[0, 'add_index'],'RI'] = np.mean(self.MGRef.at[match_table.at[0, 'add_index'],'RT_List'])
                          
                            if self.DataBase.at[i,'Tag'] =='Sample' and 'MS2_Int' in self.DataBase.at[i,'Data'].Final_Peak_Detect.keys():
                                if len(self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_Int']) > 0 and self.DataBase.at[i,'Tag']=='Sample':
                                    MS2_Int_List = self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_Int']/max(self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_Int'])
                                    MS2_MZ_List = self.DataBase.at[i, 'Data'].Final_Peak_Detect.at[ii, 'MS2_MZ']
                                    [MS2_Int_List,MS2_MZ_List] = MergeMS2(MS2_Int_List,MS2_MZ_List,self.MGRef.at[match_table.at[0, 'add_index'], 'MS2_Int'],self.MGRef.at[match_table.at[0, 'add_index'], 'MS2_MZ'])
                                    self.MGRef.at[match_table.at[0, 'add_index'], 'MS2_MZ']=MS2_MZ_List
                                    self.MGRef.at[match_table.at[0, 'add_index'], 'MS2_Int']=MS2_Int_List       
                bar.finish()
                add_RefList = pd.DataFrame({'m/z': add_MZ, 'Int': add_Int, 'RI': add_RT,self.DataBase.at[i, 'Data_Name']: add_Int, self.DataBase.at[i, 'Data_Name']+'Score': add_Score,
                                        'MS_List':add_MS_List,'RT_List':add_RT_List,'MS2_Int': add_MS2_Int, 'MS2_MZ': add_MS2_MZ})
                self.MGRef = pd.concat([self.MGRef, add_RefList])
                self.MGRef.reset_index(drop=True, inplace=True)
                self.MGRef = self.MGRef.fillna(0)
      
    def MGRawData(self,RawDataPath,RIIS_Path=''):
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
        #self.add_raw_Data(RawDataPath,Tag='Sample',Name='')
        temp_RawData = EazyMZDataProcess(RawDataPath)
        #i = len(self.DataBase)
        if RIIS_Path != '':
            #self.DataBase.at[i,'Data'].RIIS = pd.read_excel(RIIS_Path)
            temp_RawData.RIIS = pd.read_excel(RIIS_Path)
        else:
            try:
                '''
                self.DataBase.at[i,'Data'].RIIS = self.RIIS.copy()
                self.DataBase.at[i,'Data'].RIIS['RT']=0.00
                RTL = self.DataBase.at[i,'Data'].Origin_RT_List[0]
                RTR = self.DataBase.at[i,'Data'].Origin_RT_List[-1]
                for ii in range(len(self.DataBase.at[i,'Data'].RIIS)):                
                    MZ = self.DataBase.at[i,'Data'].RIIS['m/z'][ii]
                    [RT_List, Int_List]=self.DataBase.at[i,'Data'].ExtractDataPoint(MZ,RTR,RTL)
                    if max(Int_List)>50000:
                        RT_place = np.where(Int_List==max(Int_List))[0]
                        self.DataBase.at[i,'Data'].RIIS['RT'][ii] = RT_List[RT_place[0]]
                '''
                temp_RawData.RIIS = self.RIIS.copy()
                temp_RawData.RIIS['RT']=0.00
                RTL = temp_RawData.Origin_RT_List[0]
                RTR = temp_RawData.Origin_RT_List[-1]
                for ii in range(len(temp_RawData.RIIS)):                
                    MZ = temp_RawData.RIIS['m/z'][ii]
                    [RT_List, Int_List]=temp_RawData.ExtractDataPoint(MZ,RTR,RTL)
                    if max(Int_List)>50000:
                        RT_place = np.where(Int_List==max(Int_List))[0]
                        temp_RawData.RIIS['RT'][ii] = RT_List[RT_place[0]]
            except AttributeError:
                print('no RIIS in self')
                return 
            except Exception:
                print('unknown error')
                return
        temp_RawData.Final_Peak_Detect = pd.DataFrame(columns=['AverageMZ', 'RT', 'Int','RTList'])
        '''
        bar = Bar('Extrect Peak', max=len(self.MGRef))
        for ii in range(len(self.MGRef)):
            bar.next()
            MZ = self.MGRef.at[ii,'m/z']
            RI = self.MGRef.at[ii,'RI']
            RTL = RI_to_RT(RI*0.98,temp_RawData.RIIS)
            RTR = RI_to_RT(RI*1.02,temp_RawData.RIIS)
            [RT_List,Int_List] = temp_RawData.ExtractDataPoint(MZ,RTR,RTL)
            Int_List = list(map(lambda x:np.mean(Int_List[x-2:x+2]) if x in range(2,len(Int_List)-2) else Int_List[x],range(len(Int_List))))
            Peak_Place = list(filter(lambda x:Int_List[x-3]<Int_List[x-2]<Int_List[x-1]<Int_List[x] and Int_List[x+3]<Int_List[x+2]<Int_List[x+1]<Int_List[x] and Int_List[x]>7000,range(3,len(RT_List)-3)))
            if len(Peak_Place)>0:
                Peak_Place = np.where(abs(np.array(RT_List)[Peak_Place]-(RTL+RTR)/2)==min(abs(np.array(RT_List)[Peak_Place]-(RTL+RTR)/2)))[0][0]
                RT = RT_List[Peak_Place]
                Int = Int_List[Peak_Place]          
                Same_MZ_list = list(filter(lambda x:abs(temp_RawData.Final_Peak_Detect['AverageMZ'][x]-MZ)/MZ<0.000010,range(len(temp_RawData.Final_Peak_Detect))))
                Same_RT_list = list(filter(lambda x:abs(RT-temp_RawData.Final_Peak_Detect['RT'][x])<6,range(len(temp_RawData.Final_Peak_Detect))))
                Same_list = list(filter(lambda x:x in Same_MZ_list,Same_RT_list))
                if len(Same_list) > 0:
                    pass
                else:
                    temp_Final_Peak = pd.DataFrame([(MZ, RT, Int)], columns=['AverageMZ', 'RT', 'Int'])
                    temp_RawData.Final_Peak_Detect = pd.concat([temp_RawData.Final_Peak_Detect, temp_Final_Peak])
                    temp_RawData.Final_Peak_Detect.reset_index(drop=True,inplace=True)
        bar.finish()
        '''
        bar = Bar('Extrect Peak', max=len(self.MGRef))
        for ii in range(len(self.MGRef)):       
            bar.next()
            MZ = self.MGRef.at[ii,'m/z']
            RI = self.MGRef.at[ii,'RI']
            RTL = RI_to_RT(RI*0.98,temp_RawData.RIIS)
            RTR = RI_to_RT(RI*1.02,temp_RawData.RIIS)
            [RT_List,Int_List] = temp_RawData.ExtractDataPoint(MZ,RTR+6,RTL-6)
            Int_List = list(map(lambda x:np.mean(Int_List[x-2:x+2]) if x in range(2,len(Int_List)-2) else Int_List[x],range(len(Int_List))))
            Diff_List = np.diff(Int_List[1:len(Int_List)-2])
            Diff_List = np.array(list(map(lambda x:GaussSmooth(Diff_List[x-2:x+2+1]) if x in range(2,len(Diff_List)-2) else Diff_List[x],range(len(Diff_List)))))
            FD_List = np.array(list(map(lambda x: EazyMZDataProcess.TFFD(x,Int_List,RT_List), range(2, len(Int_List)-2))))
            FD_List = np.array(list(map(lambda x:GaussSmooth(FD_List[x-2:x+2+1]) if x in range(2,len(FD_List)-2) else FD_List[x],range(len(FD_List)))))
            SD_List = np.array(list(map(lambda x: EazyMZDataProcess.TFSD(x,Int_List,RT_List), range(2, len(Int_List)-2))))
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
            SD_Place = list(filter(lambda x: SD_List[x] < SD_Median, range(len(SD_List))))
            Begin_Place,Complet_BP = EazyMZDataProcess.Find_FContinuous(FD_List, FD_P, Diff_P, FDMode='P',mergeRule='Intersection',FC_Number=2)
            Begin_Place += 2
            Complet_BP += 2
            End_Place,Complet_EP = EazyMZDataProcess.Find_FContinuous(FD_List, FD_N, Diff_N, FDMode='N',mergeRule='Intersection',FC_Number=2)
            End_Place += 2
            Complet_EP += 2
            Peak_Place = EazyMZDataProcess.Find_SDChange(FD_Change, SD_Place)  # 峰顶位置
            Peak_Place += 2
            FD_Sequence = np.concatenate((Begin_Place, End_Place))
            FD_Sequence.sort()
            if len(Begin_Place) >= 1 and len(End_Place) >= 1 and Begin_Place[0] < End_Place[-1]:
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
                Peak_Top = []
                for i_top in range(len(Peak_Begin)-1, -1, -1):
                    temp_top = np.where((Peak_Begin[i_top] < Peak_Place) & (Peak_End[i_top] > Peak_Place))[0]
                    if len(temp_top) >= 1:
                        top_ran = np.array(range(len(temp_top)))
                        Int_List = np.array(Int_List)
                        top_ran = list(filter(lambda x: Int_List[Peak_Place[temp_top[x]]] == np.max(Int_List[Peak_Place[temp_top]]), top_ran))
                        Peak_Top.append(Peak_Place[temp_top[top_ran[0]]])
                    else:
                        del Peak_Begin[i_top]
                        del Peak_End[i_top]
                Peak_Top = list(map(lambda x: Peak_Top[-x], range(1, len(Peak_Top)+1)))
                for v in range(len(Peak_Begin)):
                    RT = RT_List[Peak_Top[v]]
                    Int = Int_List[Peak_Top[v]]
                    BaseValue = temp_RawData.Baseline(MZ,RT_List[Peak_End[v]+1],RT_List[Peak_Begin[v]])
                    if Int>7000 and Int/BaseValue>=3:
                        #IntList = Int_List[Peak_Begin[v]:Peak_End[v]+1]
                        RTList = RT_List[Peak_Begin[v]:Peak_End[v]+1]
                        Same_MZ_list = list(filter(lambda x:abs(temp_RawData.Final_Peak_Detect['AverageMZ'][x]-MZ)/MZ<self.__param['MS1_Tor'],range(len(self.Final_Peak_Detect))))
                        Same_RT_list = list(filter(lambda x:abs(RT-temp_RawData.Final_Peak_Detect['RT'][x])<6,range(len(temp_RawData.Final_Peak_Detect))))
                        Same_list = list(filter(lambda x:x in Same_MZ_list,Same_RT_list))
                        BaseValue = temp_RawData.Baseline(temp_RawData,MZ,RT_List[Peak_End[v]+1],RT_List[Peak_Begin[v]])
                        if len(Same_list) > 0:
                            pass
                        else:
                            temp_Final_Peak = pd.DataFrame([(MZ, RT, Int,RTList)], columns=['AverageMZ', 'RT', 'Int','RTList'])
                            temp_RawData.Final_Peak_Detect = pd.concat([temp_RawData.Final_Peak_Detect, temp_Final_Peak])
                            temp_RawData.Final_Peak_Detect.reset_index(drop=True,inplace=True)
        bar.finish()


if __name__ == '__main__':
    mp.freeze_support()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('./HeuristicPeakListGenerator.ico'))
    #main = FirstMainWindow()
    main_HPLG = HeuristicProducerUI()
    main_HPLG.show()
    sys.exit(app.exec_())
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        