#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 21:22:13 2022

@author: lacter
"""
import pymzml
import pandas as pd
from progress.bar import Bar
from matplotlib import pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


class EazyMZDataProcess(object):   
    def __init__(self,InputDataPath):
        self.file_path = InputDataPath
        run = pymzml.run.Reader(self.file_path)
        self.Time_List = []
        self.Spectrum_List = []
        for n,spectrum in enumerate(run):
            self.Time_List.append(round(spectrum.scan_time[0]/60,3))
            self.Spectrum_List.append(spectrum)
    def RTPointFind(self,RT,RT_Tol=0.05):
        ''' match RT
            example scrip：
            EazyDataProcess.RTFindRange(5.11,0.05)
            5.11 is RT,0.05 is tolerance. Means [5.05,5.16] 
            Valified_List = [RT1,RT2……]
        '''
        temp = 0
        Valified_List = []
        for i in self.Time_List:
            temp += 1
            if i > RT - RT_Tol:
                if i < RT + RT_Tol:
                    Valified_List.append(temp)
        return Valified_List
    def RTRangeFind(self,RT_Lift,RT_Right):
        ''' match RT
            example scrip：
            EazyDataProcess.RTRangeFind(5.11,6.11)
            5.11 is begin RT,6.11 is end RT, Means [5.11,6.11] 
            Valified_List = [RT1,RT2……]
        '''
        temp = 0
        Valified_List = []
        for i in self.Time_List:
            temp += 1
            if i > RT_Lift:
                if i < RT_Right:
                    Valified_List.append(temp)
        return Valified_List
    def ExportRangeData(self,
                        RT,
                        OutputPath = "RangeData.xlsx",WhetherOutput=True):
        ''' Export origin all m/z Data in selected time
            example scrip：
            EazyDataProcess.ExportRangeData(1,5)
            Data acquisition during 1-5 min, will be export.
        '''
        Valified_List = self.RTRangeFind(RT[0],RT[1])
        df = pd.DataFrame()
        bar = Bar('Processing',max=len(Valified_List))
        for i in Valified_List:
            mz_temp = []
            int_temp = []
            for ii in self.Spectrum_List[i].centroidedPeaks:
                mz_temp.append(ii[0])
                int_temp.append(ii[1])
            df_mz = pd.DataFrame({'Scan '+str(i+1):mz_temp})
            df_int = pd.DataFrame({str(self.Time_List[i]):int_temp})
            df_com = pd.concat([df_mz,df_int],axis=1)
            #df_com = df_com.drop(df_com.loc[df_com[str(self.Time_List[i])]<3000].index)
            del mz_temp
            del int_temp
            del df_mz
            del df_int
            df_temp = pd.concat([df,df_com],axis=1)
            del df_com
            del df
            df = df_temp
            del df_temp
            df.fillna(0)
            bar.next()
        bar.finish()
        if WhetherOutput != False:
            OutputPath = [str(RT[0])+"-"+str(RT[1])+OutputPath]
            df.to_excel(OutputPath[0],index=True,header=True)
    def ExportPointData_Origin(self,
                        RT,mz,
                        RT_Tol=0.1,mz_Tol=0.00001,
                        OutputPath = "-PointData.xlsx",WhetherOutput=False):
        ''' Export one selected feature with m/z and RT
            example scrip：
            EazyDataProcess.ExportPointData(426.3214,2.43)
        '''
        Valified_List = self.RTPointFind(RT,RT_Tol)
        print("ExportSelectedData")
        df = pd.DataFrame()
        bar = Bar('Processing', max=len(Valified_List))
        for i in Valified_List:
            mz_temp = []
            int_temp = []
            for ii in self.Spectrum_List[i].centroidedPeaks:
                if ii[0] > mz * (1-mz_Tol):
                    if ii[0] < mz * (1+mz_Tol):
                        mz_temp.append(ii[0])
                        int_temp.append(ii[1])
            df_mz = pd.DataFrame({'Scan '+str(i+1):mz_temp})
            df_int = pd.DataFrame({str(self.Time_List[i]):int_temp})
            df_com = pd.concat([df_mz,df_int],axis=1)
            del mz_temp
            del int_temp
            del df_mz
            del df_int
            df_temp = pd.concat([df,df_com],axis=1)
            del df_com
            del df
            df = df_temp
            del df_temp
            df.fillna(0)
            bar.next()
        bar.finish()
        if WhetherOutput != False:
            OutputPath = [str(RT)+"-"+str(mz)+OutputPath]
            df.to_excel(OutputPath[0],index=True,header=True)
        return df
    def ExportPointData(self,
                        RT,
                        mz,
                        RT_Tol=0.1,mz_Tol=0.00001,
                        OutputPath = "-PointData.xlsx",WhetherOutput=False):
        Valified_List = self.RTPointFind(RT,RT_Tol)
        print("ExportSelectedData")
        df = pd.DataFrame()
        bar = Bar('Processing', max=len(Valified_List))
        Int_List = []
        RT_List = []
        for i in Valified_List:
            Int_temp = []
            for ii in self.Spectrum_List[i].centroidedPeaks:
                if ii[0] > mz * (1-mz_Tol):
                    if ii[0] < mz * (1+mz_Tol):
                        Int_temp.append(ii[1])
            Int_List.append(max(Int_temp))
            RT_List.append(self.Time_List[i])
            bar.next() 
        bar.finish()
        Int_AverageList = []
        for i in range(2,len(Int_List)-2):
            Int_temp = (Int_List[i-2]+Int_List[i-1]+Int_List[i]+Int_List[i+1]+Int_List[i+2])/5
            Int_AverageList.append(Int_temp)
        df_Int = pd.DataFrame({str(mz):Int_AverageList})
        df_RT = pd.DataFrame({"RT":RT_List[2:len(RT_List)-2]})
        df = pd.concat([df_RT,df_Int],axis=1)
        if WhetherOutput != False:
            OutputPath = [str(mz)+OutputPath]
            df.to_excel(OutputPath[0],index=True,header=True)
        return df
    def PlotSelectedSpectrum(self):
        print("PlotSelectedSpectrum not ready for use")
        #figure = plot.figure()
        #axes = Axes3D(figure)
        x = []
        y = []
        z = []
        for i in self.Spectrum_List:    
            for ii in i.centroidedPeaks:
                x.append(ii[0])
                z.append(ii[1])
                y.append(i.scan_time)
    def ms_Slice(self,
                 mz,
                 mz_Tol = 0.00001,RT = False,WhetherOutput=False):
        ''' Extract Data of selected m/z in all RT or selected RT range
            example scrip：
            EazyDataProcess.ms_Slice(426.3214,mz_Tol = 0.000005)
            EazyDataProcess.ms_Slice(426.3214,RT = [1,5],WhetherOutput=True)
        '''
        print("ms Slice Begin")
        Int_List = []
        if RT == False:
            bar = Bar('Processing', max=len(self.Spectrum_List))
            for i in self.Spectrum_List:
                int_temp = []
                bar.next()
                for ii in i.centroidedPeaks:
                    if ii[0] > mz *(1-mz_Tol):
                        if ii[0] < mz *(1+mz_Tol):
                            int_temp.append(ii[1])
                if not int_temp:
                    Int_List.append(0)
                else:
                    Int_List.append(np.max(int_temp))
            bar.finish()
            df_RT = pd.DataFrame(self.Time_List)
            df_int = pd.DataFrame(Int_List)
            df = pd.concat([df_RT,df_int],axis=1)  
            plt.figure(figsize=(6, 3), dpi=200)
            plt.plot(self.Time_List,Int_List,0.005)
            if WhetherOutput != False:
                OutputPath = ["Slice"+str(mz)+".xlsx"]
                df.to_excel(OutputPath[0],index=False,header=False)
                plt.savefig(fname = 'fix.svg',format="svg")
        else:
            Valified_List = self.RTRangeFind(RT[0], RT[1])
            df = pd.DataFrame()
            bar = Bar('Processing', max=len(Valified_List))
            RT_temp = []
            for i in Valified_List:
                int_temp = []
                RT_temp.append(self.Time_List[i])
                for ii in self.Spectrum_List[i].centroidedPeaks:
                    if ii[0] > mz * (1-mz_Tol):
                        if ii[0] < mz * (1+mz_Tol):
                            int_temp.append(ii[1])
                if not int_temp:
                    Int_List.append(0)
                else:
                    Int_List.append(np.max(int_temp))
                bar.next()
            bar.finish()
            df_RT = pd.DataFrame(RT_temp)
            df_int = pd.DataFrame(Int_List)
            df = pd.concat([df_RT,df_int],axis=1)
            plt.figure(figsize=(6, 3), dpi=200)
            plt.plot(RT_temp,Int_List,0.005)
            if WhetherOutput != False:
                OutputPath = ["Slice"+str(mz)+".xlsx"]
                df.to_excel(OutputPath[0],index=False,header=False)
                plt.savefig(fname = 'fix.svg',format="svg")
        return df
                    
                
    
'''          
---------------------------------------------------
'''
print("Begin")
Data = EazyMZDataProcess("/Users/lacter/Documents/EazyMZDataProcess/serum-RCOOH-F-120-scan-1.mzML")
#OutputPath = "/Users/lacter/Documents/Mzml_DataProcess/Output20220701.xlsx"
#Data.ExportRangeData([3,4])    # 已验证 7.03
Data.ExportPointData(12.914,329.19669,WhetherOutput=True)  # 已验证 7.03
#a = Data.ms_Slice(329.19)
#b = Data.ms_Slice(636.12)
#print(Data.Time_List)























