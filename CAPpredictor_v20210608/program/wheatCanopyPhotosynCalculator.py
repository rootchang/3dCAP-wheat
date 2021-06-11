#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
wheatCanopyPhotosynCalculator.py

Calculate wheat canopy photosynthesis based on 3D architecture and ray tracing algorithm
"""

import argparse
import time
import numpy as np
import math
import re
import random
import copy
import os
import sys
import json
import subprocess
from datetime import datetime
#from scipy.interpolate import pchip
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit

# interpolant = PchipInterpolator(x, y)
# xnew = np.linspace(0, 10, num=41, endpoint=True)
# interpolant(xnew)

__author__="CHANG Tiangen"
__date__="20210519"
__copyright__="Copyright 2020 Aspar CHANG"
__license__="SIPPE"
__version__="0.3"
__email__="changtiangen@sippe.ac.cn"
__status__="Prototype"
__credits__=[]

class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass

def createHelp():
    """
    Create the command line interface of the programme.
    """
    
    epilog_string="Any bug is welcome reported to changtiangen@sippe.ac.cn"
    description_string='The program is going to calculate canopy photosynthesis for a given 3D wheat canopy.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-di', '--input directory', dest='dir_input', help='directory for input files')
    parser.add_argument('-i1', '--input file', dest='file_in1', help='input raw data file of manually-extracted tiller structure')
    parser.add_argument('-i2', '--input file2', dest='file_in2', type=str, help='input weather data file')
    parser.add_argument('-i3', '--input file3', dest='file_in3', type=str, nargs='+', help='input configure files (separated by space)')
    parser.add_argument('-ih', '--image height', dest='img_H', type=int, default=3264,help='image height in pixels')
    parser.add_argument('-lineNA', '--cultivar line name', dest='lineNA', type=str, help='cultivar specific key word in image file name')
    parser.add_argument('-sp', '--spike', dest='withSpike', type=int,default=1, help='Canopy with spike (1) or de-spiked (0)')
    
    parser.add_argument('-rID', '--run ID', dest='runID', type=int, help='remark of current run ID, used as file prefix and random seed.')
    parser.add_argument('-dp', '--program directory', dest='dir_prog', type=str, default='./',help='directory for programs')
    parser.add_argument('-do', '--output directory', dest='dir_result', type=str, help='directory for result output')
    parser.add_argument('-3d', '--tri', dest='file_tri', type=str, help='canopy triangulation data file')
    parser.add_argument('-rmk', '--remark', dest='remark', type=str, default='',help='Additional remark of current run, used as statistical file surffix.')
    parser.add_argument('-ow', '--overWrite', dest='ow', type=int, default=0, help='Overwrite the existing files (1) or not (0).')

    op=parser.parse_args()
    return op
    
def powerFunc(x, I0, k, Ir):
    return I0 * np.exp(-k * x) + Ir
    
def fileLineNumFastCounter(file_name):
    from itertools import (takewhile, repeat)
    buffer = 1024 * 1024
    with open(file_name) as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        return sum(buf.count('\n') for buf in buf_gen)
    
def inputFileProcessing(weatherIn):
    print('Input weather file processing...\n')
    dataIn=open(weatherIn,'r').readlines()
    DateTimeStr_list=[]
    HOUR_list=[]
    inDL_list=[]
    inSL_list=[]
    VPD_list=[] # VPD=0.61137*EXP(17.502*T/(240.97+T))*(1-RH)
    AT_list=[]
    inTL_list=[]
    
    for line_i in range(1,len(dataIn)):
        words=dataIn[line_i].strip().split('\t')
        if len(words)>10:
            DateTimeStr_list.append('-'.join(words[0:6]))
            HOUR_list.append(round(int(words[3])+int(words[4])/60+int(words[5])/3600,2))
            inDL_list.append(float(words[6]))
            inSL_list.append(float(words[7]))
            VPD_list.append(float(words[8]))
            AT_list.append(float(words[9]))
            inTL_list.append(float(words[10]))
    return DateTimeStr_list,HOUR_list,inDL_list,inSL_list,VPD_list,AT_list,inTL_list

def organPhotoCalculator(organPos,facet_area, organKt,organKr,organNC, totLight, coef_leaf_AQ_paras,coef_stem_AQ_paras,coef_grain_AQ_paras,coef_awn_AQ_paras,facetNum,AT_list):
    AT_list=np.array(AT_list)
    timePts=len(totLight[0,:])
    Amat=np.zeros([facetNum,timePts])
    organNC_old=-100
    organPos_old=-100
    facet_i=0
    while facet_i<facetNum:
        PPFD_list=totLight[facet_i,:]/(1-organKt[facet_i]-organKr[facet_i]) # convert ONLY LEAF abs to incident light
        organNC_new=organNC[facet_i]
        if organPos[facet_i]>0.5: # leaf
            coef_AQ_paras=coef_leaf_AQ_paras
        elif abs(organPos[facet_i]-0)<0.5: # stem
            coef_AQ_paras=coef_stem_AQ_paras
        elif abs(organPos[facet_i]+5)<0.5: # awn
            coef_AQ_paras=coef_awn_AQ_paras
        elif organPos[facet_i]<-1.5 and organPos[facet_i]>-4.5: # glume, branch and grain
            coef_AQ_paras=coef_grain_AQ_paras
        elif abs(organPos[facet_i]+6)<0.5: # soil
            facet_i+=1
            continue
        else:
            raise Exception('Line %d, unrecognized organPos: %d'%(facet_i+1,organPos[facet_i]))
        if abs(organNC_new-organNC_old)>0.001 or abs(organPos_old-organPos[facet_i])>0.001:
            Asat=coef_AQ_paras[0][0]*organNC_new**2+coef_AQ_paras[0][1]*organNC_new+coef_AQ_paras[0][2]
            alpha=coef_AQ_paras[1][0]*organNC_new**2+coef_AQ_paras[1][1]*organNC_new+coef_AQ_paras[1][2]
            theta=coef_AQ_paras[2][0]*organNC_new**2+coef_AQ_paras[2][1]*organNC_new+coef_AQ_paras[2][2]
            Rd=coef_AQ_paras[3][0]*organNC_new**2+coef_AQ_paras[3][1]*organNC_new+coef_AQ_paras[3][2]
            Asat_new=Asat*(coef_AQ_paras[4][0]*AT_list**2+coef_AQ_paras[4][1]*AT_list+coef_AQ_paras[4][2]) # Asat-T relation
            Rd_new=Rd*(coef_AQ_paras[5][0]*AT_list**2+coef_AQ_paras[5][1]*AT_list+coef_AQ_paras[5][2]) # Rd-T relation
            Amax=Asat_new+Rd_new
            organNC_old=organNC_new
            organPos_old=organPos[facet_i]
        A_list=(alpha*PPFD_list+Amax-np.sqrt((alpha*PPFD_list+Amax)**2-4*alpha*PPFD_list*Amax*theta))/(2*theta)-Rd_new #umol/m2/s
        #A_list=[A_list[i]*facet_area[i]/10 for i in range(timePts)] # nmol/s
        Amat[facet_i]=A_list
        facet_i+=1
    return Amat

def saveData(fnOut,header_str,data_mat,format_str,lineNum):
    if lineNum==0:
        lineNum=len(data_mat)
    
    format_str=format_str[0:-1]#.replace('\t','')
    np.savetxt(fnOut, data_mat,fmt=format_str,delimiter='\t')
    
    with open(fnOut, 'r+') as f:
        content = f.read()        
        f.seek(0, 0)
        f.write(header_str+content)
        f.close()
    
def plot_canopyPhotoStats(CP_mat,DateTimeStr_list,HOUR_list,inDL_list,inSL_list,VPD_list,AT_list,inTL_list,xmin,xmax,ymin,ymax,resultDir,fn_str,runID,fnOut_remark):
    # plot:
    #   1. Stats ear/stem/leaves/canopy/Total cumulated incident light (mol/m2 ground)
    #   2. Stats cumulated incident light on different leaf positions (mol/m2 ground)
    #   3. Scatter for average incident light of leaf facets at different canopy height (umol/m2 leaf/s)
    
    #   4. Ear/stem/leaves/canopy photosynthesis rate versus Day time (umol/m2 ground/s)
    #   5. Photosynthesis rate of different leaf positions versus Day time (umol/m2 ground/s)
    #   6. Ear/stem/leaves/canopy photosynthesis rate versus total ambient incident PPFD (umol/m2 ground/s)
    #   7. Photosynthesis rate of different leaf positions versus total ambient incident PPFD (umol/m2 ground/s)
    #   8. Scatter for average photosynthesis rate of leaf facets at different canopy height (umol/m2 leaf/s)
    
    #   9. Stats ear/stem/leaves/canopy cumulated photosynthesis (mol/m2 ground)
    #   10. Stats photosynthesis on different leaf positions (mol/m2 ground)
    
    # fonts setting
    ticklabelFont=10
    labelFont=12
    legendFont=10
    refFont=8
    
    x_time = [datetime.strptime(c, "%Y-%m-%d-%H-%M-%S") for c in DateTimeStr_list]
    
    tot_timePts=len(HOUR_list)
    timeInterval=abs(np.diff(HOUR_list))
    timeInterval=np.hstack(([0],timeInterval,[0]))
    timePeriod_list=[(timeInterval[i+1]+timeInterval[i])/2 for i in range(len(timeInterval)-1)] # for time integration
    tot_facetNum=len(CP_mat)
    facetOrganCode=CP_mat[:,3]
    facetOrganPos=CP_mat[:,2]
    facetArea=CP_mat[:,17]
    facetArea_mat=np.tile(CP_mat[:,17:18],[1,tot_timePts]) # cm2
    facetPPFD=CP_mat[:,17+7:17+7*(tot_timePts+1):7] # incident PPFD on facet (umol/m2/s)
    facetA=CP_mat[:,17+7*tot_timePts+1:] # A of facet (umol/m2/s)
    facetPPF=facetArea_mat*facetPPFD/10 # nmol/s
    facetCA=facetArea_mat*facetA/10 # nmol/s
    
    stem_facet_logic=(facetOrganCode==3)+0
    leaf_facet_logic=(facetOrganCode==4)+0
    ear_facet_logic=(facetOrganCode==5)+0
    awn_facet_logic=(facetOrganCode==5)*(facetOrganPos==-5)
    nonawn_facet_logic=(facetOrganCode==5)*(facetOrganPos>-5)
    
    # calculate leaf facet PPFD distribution at different time points
    leaf_facet_area=facetArea*leaf_facet_logic/((xmax-xmin)*(ymax-ymin))
    PPFD_binWid=100
    PPFD_upperBound=1100
    PPFD_grad=np.array(range(0,PPFD_upperBound,PPFD_binWid))
    facet_area_PPFD=[[0]*len(PPFD_grad) for _ in range(tot_timePts)]
    for time_i in range(tot_timePts):
        facetPPFD_temp=CP_mat[:,17+7*(time_i+1)]
        for PPFD_j in range(1,len(PPFD_grad)):
            facet_temp=((facetPPFD_temp<PPFD_grad[PPFD_j])+0)*((facetPPFD_temp>=PPFD_grad[PPFD_j-1])+0)
            facet_area_PPFD[time_i][PPFD_j-1]=sum(leaf_facet_area*facet_temp) # m2/m2
        facet_temp=((facetPPFD_temp>=PPFD_grad[PPFD_j])+0)
        facet_area_PPFD[time_i][PPFD_j]=sum(leaf_facet_area*facet_temp) # m2/m2
        
    max_pos=int(max(facetOrganPos))
    leaf_i_facet_logic=np.zeros([max_pos,tot_facetNum])
    for i in range(1,max_pos+1):
        leaf_i_facet_logic[i-1]=(facetOrganCode==4)*(facetOrganPos==i)+0
    
    tot_leaf_facet_area=sum(leaf_facet_area)
    tot_stem_facet_area=sum(facetArea*stem_facet_logic)/((xmax-xmin)*(ymax-ymin))
    tot_ear_facet_area=sum(facetArea*ear_facet_logic)/((xmax-xmin)*(ymax-ymin))
    print('LAI, CAI, EAI: ',tot_leaf_facet_area,tot_stem_facet_area,tot_ear_facet_area) # m2/m2 ground
    
    leaf_timePtsPPF=np.dot(leaf_facet_logic,facetPPF)/1000/((xmax-xmin)*(ymax-ymin)/10000) # all-leaf PPF at different timePts, umol/m2 ground/s
    stem_timePtsPPF=np.dot(stem_facet_logic,facetPPF)/1000/((xmax-xmin)*(ymax-ymin)/10000) # all-stem PPF at different timePts, umol/m2 ground/s
    ear_timePtsPPF=np.dot(ear_facet_logic,facetPPF)/1000/((xmax-xmin)*(ymax-ymin)/10000) # all-ear PPF at different timePts, umol/m2 ground/s
    awn_timePtsPPF=np.dot(awn_facet_logic,facetPPF)/1000/((xmax-xmin)*(ymax-ymin)/10000) # awn PPF at different timePts, umol/m2 ground/s
    nonawn_timePtsPPF=np.dot(nonawn_facet_logic,facetPPF)/1000/((xmax-xmin)*(ymax-ymin)/10000) # non-awn PPF at different timePts, umol/m2 ground/s
    canopy_timePtsPPF=list(leaf_timePtsPPF+stem_timePtsPPF+ear_timePtsPPF) # umol/m2 ground/s at different time points
    
    leaf_timePtsCA=np.dot(leaf_facet_logic,facetCA)/1000/((xmax-xmin)*(ymax-ymin)/10000) # all-leaf CA at different timePts, umol/m2 ground/s
    stem_timePtsCA=np.dot(stem_facet_logic,facetCA)/1000/((xmax-xmin)*(ymax-ymin)/10000) # all-stem CA at different timePts, umol/m2 ground/s
    ear_timePtsCA=np.dot(ear_facet_logic,facetCA)/1000/((xmax-xmin)*(ymax-ymin)/10000) # all-ear CA at different timePts, umol/m2 ground/s
    awn_timePtsCA=np.dot(awn_facet_logic,facetCA)/1000/((xmax-xmin)*(ymax-ymin)/10000) # awn CA at different timePts, umol/m2 ground/s
    nonawn_timePtsCA=np.dot(nonawn_facet_logic,facetCA)/1000/((xmax-xmin)*(ymax-ymin)/10000) # non-awn CA at different timePts, umol/m2 ground/s
    canopy_timePtsCA=list(leaf_timePtsCA+stem_timePtsCA+ear_timePtsCA) # umol/m2 ground/s at different time points
    
    leaf_i_timePtsPPF=np.zeros([max_pos,tot_timePts])
    leaf_i_timePtsCA=np.zeros([max_pos,tot_timePts])
    for i in range(1,max_pos+1):
        leaf_i_timePtsPPF[i-1]=np.dot(leaf_i_facet_logic[i-1],facetPPF)/1000/((xmax-xmin)*(ymax-ymin)/10000) # leaf_i PPF at different timePts, umol/m2 ground/s
        leaf_i_timePtsCA[i-1]=np.dot(leaf_i_facet_logic[i-1],facetCA)/1000/((xmax-xmin)*(ymax-ymin)/10000) # leaf_i CA at different timePts, umol/m2 ground/s
    
    # PPFD versus canopy height
    layer_thickness=1 # cm
    max_canopy_height=CP_mat[:,[7,10,13]].max()
    canopy_layers=int(max_canopy_height//layer_thickness+1)
    canopy_depth_array_temp=[layer_thickness]*canopy_layers
    canopy_depth_array=[sum(canopy_depth_array_temp[0:i]) for i in range(canopy_layers)]
#     canopy_height_list=copy.deepcopy(canopy_depth_array)
#     canopy_height_list.reverse()
#     canopy_height_list=np.array(canopy_height_list)
    canopy_depth_array=np.array(canopy_depth_array)
    canopy_GAI_array=np.zeros(canopy_layers)
    canopy_PPF_abs_mat=np.zeros([canopy_layers,tot_timePts])
    for i in range(tot_facetNum):
        layer_temp=int(np.mean(CP_mat[i,[7,10,13]])//layer_thickness+1)
        PPFD_temp=facetPPF[i,:]/((xmax-xmin)*(ymax-ymin)/10000) # nmol/m2 ground/s
        canopy_PPF_abs_mat[canopy_layers-layer_temp,:]+=PPFD_temp
        canopy_GAI_array[canopy_layers-layer_temp]+=facetArea[i] # cm2
    canopy_PPF_abs_mat=canopy_PPF_abs_mat/1000 # umol/m2 ground/s
    canopy_PPF_abs_mat_cum=np.zeros([canopy_layers,tot_timePts])
    for i in range(tot_timePts):
        canopy_PPF_abs_temp=canopy_PPF_abs_mat[:,i]
        canopy_PPF_abs_mat_cum[:,i]=np.cumsum(canopy_PPF_abs_temp)
    canopy_PPF_remain_mat=np.tile(inTL_list,[canopy_layers,1]) - canopy_PPF_abs_mat_cum # umol/m2 ground/s
    canopy_GAI_array=canopy_GAI_array/((xmax-xmin)*(ymax-ymin)) #m2/m2 ground
    canopy_GAI_array_cum=np.cumsum(canopy_GAI_array)
    
    
    # CP versus Time
    plt.figure(figsize=(5.2,3))
    left_edge=0.2
    bottom_edge=0.2
    width=1-left_edge-0.02
    height=1-bottom_edge-0.02
    wid_space=0.05
    hei_space=0.05
    row_num=1
    col_num=1
    panel_width=(width-wid_space*(col_num-1))/col_num
    panel_height=(height-hei_space*(row_num-1))/row_num
    rect_line=[]
    ax_list=[]
    for i in range(row_num):
        for j in range(col_num):
            rect_line.append([left_edge+j*(panel_width+wid_space),bottom_edge+(row_num-i-1)*(panel_height+hei_space),\
                              panel_width,panel_height])
    for i in range(len(rect_line)):
        ax_list.append(plt.axes(rect_line[i]))
        ax_list[i].tick_params(direction='in',axis='both',length=2, which='major', labelsize=ticklabelFont,pad=2,top='off', right='off')
        ax_list[i].yaxis.set_label_coords(-0.1, 0.5)
        ax_list[i].xaxis.set_label_coords(0.5, -0.15)
    ax_list[0].plot(x_time,canopy_timePtsCA,linestyle='-',linewidth=2,color='k',label='Canopy')
    ax_list[0].plot(x_time,leaf_timePtsCA,linestyle='-',linewidth=1,color='g',label='Leaf')
    ax_list[0].plot(x_time,stem_timePtsCA,linestyle='-',linewidth=1,color='grey',label='Culm')
    ax_list[0].plot(x_time,ear_timePtsCA,linestyle='-',linewidth=1,color='y',label='Ear')
    ax_list[0].legend(loc=(0.0,0.7),fontsize=legendFont,frameon=False,labelspacing=0,handlelength=1)
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gcf().autofmt_xdate()
    ax_list[0].set_xlabel('Daytime')
    ax_list[0].set_ylabel(r'$\mathregular{A_{c\ net}\ (\mu mol\ m^{-2}\ ground\ s^{-1})}$')
    #plt.show()
    fig_file_name=os.path.join(resultDir,'Run'+str(runID)+'_CP_PPFD_'+str(fn_str)+'_'+fnOut_remark+'.png')
    plt.savefig(fig_file_name,format="png",dpi=600)
    plt.close()
    
    # PPFD versus canopy height
    plt.figure(figsize=(5.2,3))
    left_edge=0.2
    bottom_edge=0.2
    width=1-left_edge-0.02
    height=1-bottom_edge-0.02
    wid_space=0.05
    hei_space=0.05
    row_num=1
    col_num=1
    panel_width=(width-wid_space*(col_num-1))/col_num
    panel_height=(height-hei_space*(row_num-1))/row_num
    rect_line=[]
    ax_list=[]
    for i in range(row_num):
        for j in range(col_num):
            rect_line.append([left_edge+j*(panel_width+wid_space),bottom_edge+(row_num-i-1)*(panel_height+hei_space),\
                              panel_width,panel_height])
    for i in range(len(rect_line)):
        ax_list.append(plt.axes(rect_line[i]))
        ax_list[i].tick_params(direction='in',axis='both',length=2, which='major', labelsize=ticklabelFont,pad=2,top='off', right='off')
        ax_list[i].yaxis.set_label_coords(-0.1, 0.5)
        ax_list[i].xaxis.set_label_coords(0.5, -0.15)
    ax_list[0].plot(x_time,canopy_timePtsCA,linestyle='-',linewidth=2,color='k',label='Canopy')
    ax_list[0].plot(x_time,leaf_timePtsCA,linestyle='-',linewidth=1,color='g',label='Leaf')
    ax_list[0].plot(x_time,stem_timePtsCA,linestyle='-',linewidth=1,color='grey',label='Culm')
    ax_list[0].plot(x_time,ear_timePtsCA,linestyle='-',linewidth=1,color='y',label='Ear')
    ax_list[0].legend(loc=(0.0,0.7),fontsize=legendFont,frameon=False,labelspacing=0,handlelength=1)
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gcf().autofmt_xdate()
    ax_list[0].set_xlabel('Daytime')
    ax_list[0].set_ylabel(r'$\mathregular{A_{c\ net}\ (\mu mol\ m^{-2}\ ground\ s^{-1})}$')
    #plt.show()
    fig_file_name=os.path.join(resultDir,'Run'+str(runID)+'_CP_PPFD_'+str(fn_str)+'_'+fnOut_remark+'.png')
    plt.savefig(fig_file_name,format="png",dpi=600)
    plt.close()
    
    
    # Total area in different PPFD intervals barplot
    plt.figure(figsize=(5.2,3))
    left_edge=0.2
    bottom_edge=0.2
    width=1-left_edge-0.02
    height=1-bottom_edge-0.02
    wid_space=0.05
    hei_space=0.05
    row_num=1
    col_num=1
    panel_width=(width-wid_space*(col_num-1))/col_num
    panel_height=(height-hei_space*(row_num-1))/row_num
    rect_line=[]
    ax_list=[]
    for i in range(row_num):
        for j in range(col_num):
            rect_line.append([left_edge+j*(panel_width+wid_space),bottom_edge+(row_num-i-1)*(panel_height+hei_space),\
                              panel_width,panel_height])
    for i in range(len(rect_line)):
        ax_list.append(plt.axes(rect_line[i]))
        ax_list[i].tick_params(direction='in',axis='both',length=2, which='major', labelsize=ticklabelFont,pad=2,top='off', right='off')
        ax_list[i].yaxis.set_label_coords(-0.1, 0.5)
        ax_list[i].xaxis.set_label_coords(0.5, -0.15)
    width=100 # umol/m2/s
    index_query=HOUR_list.index(12)
    timePt_query=HOUR_list[index_query]
    ear_timePtsPPF_plot=ear_timePtsPPF[index_query]
    ax_list[0].bar(PPFD_grad+PPFD_binWid/2, facet_area_PPFD[index_query], width=width,fc = 'k')
    ax_list[0].scatter([ear_timePtsPPF_plot],[0.5],c = 'y',marker = 'o')
    ax_list[0].set_xlabel(r'$\mathregular{PPFD\ (\mu mol\ m^{-2}\ ground\ s^{-1})}$')
    ax_list[0].set_ylabel(r'$\mathregular{Leaf\ area\ index\ (m^2\ m^{-2})}$')
    ax_list[0].set_xlim([0,PPFD_upperBound+PPFD_binWid/2])
    ax_list[0].set_ylim([0,round(facet_area_PPFD[index_query][0]*1.2)])
    #plt.show()
    fig_file_name=os.path.join(resultDir,'Run'+str(runID)+'_Area_PPFD_'+str(fn_str)+'_'+str(timePt_query)+'_'+fnOut_remark+'.png')
    plt.savefig(fig_file_name,format="png",dpi=600)
    plt.close()
    # save data to txt
    fig_file_name=os.path.join(resultDir,'Run'+str(runID)+'_PPFD_LAI_'+str(fn_str)+'_'+fnOut_remark+'.txt')
    fnOut=open(fig_file_name,'w')
    fnOut.write('Time(oclock)/PPFD(umol/m2/s)'+'\t'+'\t'.join([str(PPFD_grad[i_temp])+'-'+str(PPFD_grad[i_temp+1]) for i_temp in range(len(PPFD_grad)-1)])+'\t>'+str(PPFD_grad[-1])+'\n')
    for index_query in range(len(HOUR_list)):
        timePt_query=HOUR_list[index_query]
        fnOut.write(str(timePt_query)+'\t'+'\t'.join([str(c) for c in facet_area_PPFD[index_query]])+'\n')
    fnOut.close()
    
    # Canopy radiation extinction profile plot
    plt.figure(figsize=(6.4,3))
    left_edge=0.1
    bottom_edge=0.2
    width=1-left_edge-0.02
    height=1-bottom_edge-0.02
    wid_space=0.1
    hei_space=0.05
    row_num=1
    col_num=2
    panel_width=(width-wid_space*(col_num-1))/col_num
    panel_height=(height-hei_space*(row_num-1))/row_num
    rect_line=[]
    ax_list=[]
    for i in range(row_num):
        for j in range(col_num):
            rect_line.append([left_edge+j*(panel_width+wid_space),bottom_edge+(row_num-i-1)*(panel_height+hei_space),\
                              panel_width,panel_height])
    for i in range(len(rect_line)):
        ax_list.append(plt.axes(rect_line[i]))
        ax_list[i].tick_params(direction='in',axis='both',length=2, which='major', labelsize=ticklabelFont,pad=2,top='off', right='off')
        ax_list[i].yaxis.set_label_coords(-0.15, 0.5)
        ax_list[i].xaxis.set_label_coords(0.5, -0.1)
    index_query=HOUR_list.index(12)
    timePt_query=HOUR_list[index_query]
    canopy_PPF_remain=canopy_PPF_remain_mat[:,index_query]
    #ax_list[0].plot(canopy_GAI_array_cum,canopy_PPF_remain,c='k', ls='-', lw=2) # ,marker='o', mec='k', mfc='k',ms=10
    ax_list[0].scatter(canopy_GAI_array_cum,canopy_PPF_remain,c = 'k',marker = 'o',s=5)
    ax_list[0].set_xlabel(r'$\mathregular{Cummulative\ GAI}$')
    ax_list[0].set_ylabel(r'$\mathregular{PPFD_{remained}\ (\mu mol\ m^{-2}\ s^{-1})}$')
    PPFD_incident=inTL_list[index_query]
    ax_list[0].set_xlim([0,round(canopy_GAI_array_cum[-1]*1.2)])
    ax_list[0].set_ylim([0,round(0.01*PPFD_incident)*100*1.2])
    
    #I0=canopy_PPF_remain[0]-canopy_PPF_remain[-1]
    #Ir=canopy_PPF_remain[-1]
    #custom_powerFunc = lambda x, k: powerFunc(x, I0, k, Ir)
#     print(canopy_GAI_array_cum,canopy_PPF_remain)
#     fnOut_temp=open('test.txt','w')
#     fnOut_temp.write('\t'.join([str(c) for c in canopy_GAI_array_cum])+'\n')
#     fnOut_temp.write('\t'.join([str(c) for c in canopy_PPF_remain])+'\n')
#     fnOut_temp.close()
    
    popt, pcov = curve_fit(powerFunc, canopy_GAI_array_cum, canopy_PPF_remain,bounds=np.array([(0,0,0),(10000,1,10000)])) #in list popt, it is the extinction coef. k
    y2 = [powerFunc(i, popt[0], popt[1], popt[2]) for i in canopy_GAI_array_cum]
    ax_list[0].plot(canopy_GAI_array_cum,y2,'r-',lw=1)
    if popt[2]<0.001:
        ax_list[0].text(0.5,canopy_PPF_remain[0]*1.1,r'$\mathregular{I\ =\ %d \times e^{-%.2f\cdot LAI}}$'%(popt[0],popt[1]),ha = 'left',va = 'top',fontsize=legendFont)
    else:
        ax_list[0].text(0.5,canopy_PPF_remain[0]*1.1,r'$\mathregular{I\ =\ %d \times e^{-%.2f\cdot LAI}+%d}$'%(popt[0],popt[1],popt[2]),ha = 'left',va = 'top',fontsize=legendFont)
    ax_list[0].text(0.5,canopy_PPF_remain[0]*0.3,r'$\mathregular{LAI\ =\ %.1f}$'%tot_leaf_facet_area,ha = 'left',va = 'top',fontsize=refFont)
    ax_list[0].text(0.5,canopy_PPF_remain[0]*0.2,r'$\mathregular{CAI\ =\ %.1f}$'%tot_stem_facet_area,ha = 'left',va = 'top',fontsize=refFont)
    ax_list[0].text(0.5,canopy_PPF_remain[0]*0.1,r'$\mathregular{EAI\ =\ %.1f}$'%tot_ear_facet_area,ha = 'left',va = 'top',fontsize=refFont)
    #ax_list[1].plot(canopy_depth_array - layer_thickness/2,canopy_PPF_remain,c='k', ls='-', lw=2) # ,marker='o', mec='k', mfc='k',ms=10
    ax_list[1].scatter(canopy_depth_array - layer_thickness/2,canopy_PPF_remain,c = 'k',marker = 'o',s=5)
    ax_list[1].set_xlabel(r'$\mathregular{Canopy\ depth\ (cm)}$')
    ax_list[1].set_ylabel(r'$\mathregular{PPFD_{remained}\ (\mu mol\ m^{-2}\ s^{-1})}$')
    PPFD_incident=inTL_list[index_query]
    ax_list[1].set_xlim([0,round(0.1*canopy_depth_array[-1])*10*1.2])
    ax_list[1].set_ylim([0,round(0.01*PPFD_incident)*100*1.2])
    
    #plt.show()
    fig_file_name=os.path.join(resultDir,'Run'+str(runID)+'_Height_absPPFD_'+str(fn_str)+'_'+str(timePt_query)+'_'+fnOut_remark+'.png')
    plt.savefig(fig_file_name,format="png",dpi=600)
    plt.close()
        
    # save data to txt
    fig_file_name=os.path.join(resultDir,'Run'+str(runID)+'_PPFD_depth_'+str(fn_str)+'_'+fnOut_remark+'.txt')
    fnOut=open(fig_file_name,'w')
    fnOut.write('Time(oclock)/depth(cm)'+'\t'+'\t'.join([str(c) for c in canopy_depth_array - layer_thickness/2])+'\n')
    for index_query in range(len(HOUR_list)):
        timePt_query=HOUR_list[index_query]
        canopy_PPF_remain=canopy_PPF_remain_mat[:,index_query]
        fnOut.write(str(timePt_query)+'\t'+'\t'.join([str(c) for c in canopy_PPF_remain])+'\n')
    fnOut.close()
    
    CP_PPFD_fn=os.path.join(resultDir,'Run'+str(runID)+'_CP_PPFD_'+str(fn_str)+'_'+fnOut_remark+'.txt')
    header_str_CP_PPFD='PPFD_incident(umol/m2/s)\tTime(oclock)\tPPFD_canopy\tAcanopy_net\tPPFD_leaf\tAleaf_net\tPPFD_ear\tAear_net\tPPFD_awn\tAawn_net\tPPFD_nonawn\tAnonawn_net\tPPFD_stem\tAstem_net\t'
    CP_PPFD_mat=[inTL_list,HOUR_list,canopy_timePtsPPF,canopy_timePtsCA,leaf_timePtsPPF,leaf_timePtsCA,ear_timePtsPPF,ear_timePtsCA,awn_timePtsPPF,awn_timePtsCA,nonawn_timePtsPPF,nonawn_timePtsCA,stem_timePtsPPF,stem_timePtsCA]
    for i in range(0,max_pos):
        CP_PPFD_mat=CP_PPFD_mat+[leaf_i_timePtsPPF[i]]+[leaf_i_timePtsCA[i]]
        header_str_CP_PPFD+='PPFD_leaf'+str(i+1)+'\tAleaf'+str(i+1)+'_net\t'
    header_str_CP_PPFD=header_str_CP_PPFD[0:-1]+'\n'
    col_num=len(CP_PPFD_mat)
    CP_PPFD_mat=np.array(CP_PPFD_mat).T
    format_CP_PPFD='%.2f\t'*col_num
    saveData(CP_PPFD_fn,header_str_CP_PPFD,CP_PPFD_mat,format_CP_PPFD,0)
    
def canopyReconstructor(canopyConstructorLoc,inputDir,keyPoints_raw_fn_short,M_plant_fn,M_plant_clean_fn,img_H,cultivarID,runID,config_morphology_fn,config_physiology_fn,config_rayTracing_fn,overWrite,withSpike):
    ### read in config file of canopy structure reconstruction
    # Parameters initiation
    for line in open(config_rayTracing_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='xmin':
            xmin=float(words[1])
        elif words[0]=='xmax':
            xmax=float(words[1])
        elif words[0]=='ymin':
            ymin=float(words[1])
        elif words[0]=='ymax':
            ymax=float(words[1])
        elif words[0]=='zmin':
            zmin=float(words[1])
        elif words[0]=='zmax':
            zmax=float(words[1])
    config_str=config_morphology_fn.split('\\')[-1] + ' ' + config_physiology_fn.split('\\')[-1]
    command='python %s -di %s -i %s -i2 %s -o %s -ih %d -lineNA %s -rID %d -log 1 -ow %d -sp %d'%(canopyConstructorLoc,inputDir,keyPoints_raw_fn_short,config_str,M_plant_fn,img_H,cultivarID,runID,overWrite,withSpike)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)
    print('Canopy ray tracing plant structural file preparing ...\n')
    #### read in config file of ray tracing
    M_plant_raw=open(M_plant_fn,'r')
    M_plant_clean_fh=open(M_plant_clean_fn,'w')
    for line in M_plant_raw:
        words=[float(c) for c in line.strip().split('\t')]
        judge1=(words[5]<xmax and words[5]>xmin and words[8]<xmax and words[8]>xmin and words[11]<xmax and words[11]>xmin) #x
        judge2=(words[6]<ymax and words[6]>ymin and words[9]<ymax and words[9]>ymin and words[12]<ymax and words[12]>ymin) #x
        judge3=(words[7]<zmax and words[7]>zmin and words[10]<zmax and words[10]>zmin and words[13]<zmax and words[13]>zmin) #x
        if judge1 and judge2 and judge3:
            M_plant_clean_fh.write(line)
    M_plant_clean_fh.close()

def fieldRayTracer(M_plant_clean_fn,config_rayTracing_fn,weatherfnOutDL,weatherfnOutSL,fnOutDL,fnOutSL,fastTracerLoc,current_DOY,DL_timePts,format_DLfn,format_SLfn):
    # Parameters initiation
    for line in open(config_rayTracing_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='latitude':
            lat_field=float(words[1])
        elif words[0]=='start_time':
            start_time=int(words[1])
        elif words[0]=='end_time':
            end_time=int(words[1])
        elif words[0]=='time_interv':
            time_interv=int(words[1])
        elif words[0]=='xmin':
            xmin=float(words[1])
        elif words[0]=='xmax':
            xmax=float(words[1])
        elif words[0]=='ymin':
            ymin=float(words[1])
        elif words[0]=='ymax':
            ymax=float(words[1])
        elif words[0]=='zmin':
            zmin=float(words[1])
        elif words[0]=='zmax':
            zmax=float(words[1])
        elif words[0]=='rayDense':
            rayDense=float(words[1])
    print('Canopy ray tracing for direct light ...\n')
    command='%s -L %.2f -d %d -D %.1f %.1f %.1f %.1f %.1f %.1f -S 12 -n %.2f -W %.3f %.3f %.3f -m %s -o %s -C %s -s 1'%(fastTracerLoc,lat_field,current_DOY,xmin,xmax,ymin,ymax,zmin,zmax,rayDense,start_time,time_interv,end_time,M_plant_clean_fn,fnOutDL,weatherfnOutDL)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)
    print('Canopy ray tracing for scatter light ...\n')
    command='%s -L %.2f -d %d -D %.1f %.1f %.1f %.1f %.1f %.1f -S 12 -n %.2f -W 6 1 6 -m %s -o %s -C %s -s 1'%(fastTracerLoc,lat_field,current_DOY,xmin,xmax,ymin,ymax,zmin,zmax,rayDense,M_plant_clean_fn,fnOutSL,weatherfnOutSL)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)

def DLSL_extractor(fnOutDL,fnOutSL,weatherFile,config_rayTracing_fn):
    # Parameters initiation
    for line in open(config_rayTracing_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='start_time':
            start_time=int(words[1])
        elif words[0]=='end_time':
            end_time=int(words[1])
        elif words[0]=='time_interv':
            time_interv=int(words[1])
    DateTimeStr_list,HOUR_list,inDL_list,inSL_list,VPD_list,AT_list,inTL_list=inputFileProcessing(weatherFile)
    given_timePts_num=len(HOUR_list)
    DL_timePts=int((end_time-start_time)/time_interv+1)
    
    print('Loading and processing standard DL/SL ray tracing results! the DL part ...\n')
    facetNum=fileLineNumFastCounter(fnOutDL)-1 # has header line
    colNum_fnOutDL=5+9+4+DL_timePts*7
    alldataDL=np.zeros([facetNum,colNum_fnOutDL])
    line_i=0
    with open(fnOutDL,'r') as f:
        f.readline()
        for line in f:
            line=line.replace('-nan(ind)','0')
            content_temp=[float(c) for c in line.strip().split('\t')]
            alldataDL[line_i]=content_temp
            line_i+=1
    canopy_mat=alldataDL[:,0:18]
    DL_onOrgan=alldataDL[:,-DL_timePts*7:]
    DL_time=list(range(start_time,end_time+1,time_interv))
    DL_onOrgan_time=np.zeros([facetNum,given_timePts_num*7])
    for time_i in range(given_timePts_num):
        DL_temp=inDL_list[time_i]
        if inTL_list[time_i]<10: # total light, umol/m2/s
            continue
        if HOUR_list[time_i]<start_time or HOUR_list[time_i]>end_time:
            time_temp_now=24-HOUR_list[time_i]
            if time_temp_now<start_time or time_temp_now>end_time:
                raise Exception('Photosynthsis at %s in weather data file cannot be simulated! Adjusting start_time and end_time in config_rayTracing.txt!'%DateTimeStr_list[time_i])
        else:
            time_temp_now=HOUR_list[time_i]
        if time_temp_now in DL_time:
            DL_time_index_temp= DL_time.index(time_temp_now)
            DL_onOrgan_now=DL_onOrgan[:,DL_time_index_temp*7:(DL_time_index_temp+1)*7]*DL_temp/1000
        else:
            left_time_index=int(np.floor((time_temp_now-start_time)/time_interv))
            right_time_index=int(np.ceil((time_temp_now-start_time)/time_interv))
            left_dist=np.mod(time_temp_now-start_time,time_interv)/time_interv
            DL_onOrgan_now=(DL_onOrgan[:,left_time_index*7:(left_time_index+1)*7]*(1-left_dist)+DL_onOrgan[:,right_time_index*7:(right_time_index+1)*7]*left_dist)*DL_temp/1000
        DL_onOrgan_time[:,time_i*7:(time_i+1)*7]=DL_onOrgan_now
    
    print('Loading and processing standard DL/SL ray tracing results! the SL part ...\n')
    SL_onOrgan_time=np.zeros([facetNum,given_timePts_num*7])
    inSL_array=np.array(inSL_list).reshape(given_timePts_num,1)
    colNum_fnOutSL=5+9+4+7
    line_i=0
    with open(fnOutSL,'r') as f:
        f.readline()
        for line in f:
            line=line.replace('-nan(ind)','0')
            content_temp=[float(c) for c in line.strip().split('\t')]
            dataSL=np.array(content_temp[-7:])
            SL_onOrgan_temp=np.dot(inSL_array,dataSL.reshape(1,7))/1000
            SL_onOrgan_time[line_i,:]=SL_onOrgan_temp.reshape(1,7*given_timePts_num)
            line_i+=1
    return DL_onOrgan_time,SL_onOrgan_time,canopy_mat,facetNum

def fn_preprocessing(fn_str,runID):
    fn_list=fn_str.split('\\')
    resultDir='\\'.join(fn_list[0:-1])
    if not resultDir:
        resultDir='.'
    if not os.path.exists(resultDir):
        print('  WARNING: Directory %s not exists, creating it ...'%fn_str)
        try:
            os.mkdir(resultDir)
        except:
            raise Exception('ERROR: Cannot create directory: %s'%resultDir)
    if 'Run' not in fn_list[-1]:
        fn_list[-1]='Run'+str(runID)+'_'+fn_list[-1]
    fnOut='\\'.join(fn_list)
    return fnOut,fn_list[-1]

if __name__=="__main__":
    start=time.time()
    #path = os.path.abspath(os.path.dirname(__file__))
    #type = sys.getfilesystemencoding()
    #sys.stdout = Logger('./log_wheatCanopyPhotosynCalculator')
    #command_line_str='python '
    width_format=100
    welcome_info='Welcome to use wheat canopy photosynthesis analyzer'
    time_info=time.asctime( time.localtime(time.time()) )
    print(welcome_info.center(width_format, '-'))
    print(time_info.center(width_format, ' '))
    
    op=createHelp()
    runID=op.runID
    random.seed(runID)
    inputDir=op.dir_input
    progDir=op.dir_prog 
    resultDir=op.dir_result
    keyPoints_raw_fn_short=op.file_in1
    M_plant_fn_fake_short=op.file_tri
    M_plant_fn_fake=os.path.join(resultDir,M_plant_fn_fake_short)
    M_plant_fn,M_plant_fn_short=fn_preprocessing(M_plant_fn_fake,runID)
    op.file_tri=M_plant_fn
    M_plant_clean_fn=M_plant_fn[0:-4]+'_rayTracing'+M_plant_fn[-4:]
    weatherFile=os.path.join(inputDir,op.file_in2)
    config_morphology_fn,config_physiology_fn,config_rayTracing_fn,config_canopyPhotosynthesis_fn=op.file_in3
    config_morphology_fn=os.path.join(inputDir,config_morphology_fn)
    config_physiology_fn=os.path.join(inputDir,config_physiology_fn)
    config_rayTracing_fn=os.path.join(inputDir,config_rayTracing_fn)
    config_canopyPhotosynthesis_fn=os.path.join(inputDir,config_canopyPhotosynthesis_fn)
    img_H=op.img_H
    cultivarID=op.lineNA
    fnOut_remark=op.remark
    overWrite=op.ow
    withSpike=op.withSpike
    fastTracerLoc=os.path.join(progDir,'fastTracerVS2019.exe')
    canopyConstructorLoc=os.path.join(progDir,'wheatCanopyConstructor2.py')
    temp_dir=os.path.join(resultDir,'fastTracerTemp')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    
    data = json.dumps(vars(op), indent=4,ensure_ascii=False, sort_keys=False,separators=(',', ':'))
    print('Running parameters:')
    print(data)
    
    coef_leaf_Asat_coef=1
    coef_leaf_alpha_coef=1
    coef_leaf_theta_coef=1
    coef_leaf_Rd_coef=1
    coef_culm_Asat_coef=1
    coef_culm_alpha_coef=1
    coef_culm_theta_coef=1
    coef_culm_Rd_coef=1
    coef_grain_Asat_coef=1
    coef_grain_alpha_coef=1
    coef_grain_theta_coef=1
    coef_grain_Rd_coef=1
    coef_awn_Asat_coef=1
    coef_awn_alpha_coef=1
    coef_awn_theta_coef=1
    coef_awn_Rd_coef=1
    # Parameters initiation
    for line in open(config_physiology_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='coef_leaf_Asat':
            coef_leaf_Asat=[float(c) for c in words[1:]]
        elif words[0]=='coef_leaf_alpha':
            coef_leaf_alpha=[float(c) for c in words[1:]]
        elif words[0]=='coef_leaf_theta':
            coef_leaf_theta=[float(c) for c in words[1:]]
        elif words[0]=='coef_leaf_Rd':
            coef_leaf_Rd=[float(c) for c in words[1:]]
        elif words[0]=='coef_leaf_Asat_airT':
            coef_leaf_Asat_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_leaf_Rd_airT':
            coef_leaf_Rd_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_culm_Asat':
            coef_culm_Asat=[float(c) for c in words[1:]]
        elif words[0]=='coef_culm_alpha':
            coef_culm_alpha=[float(c) for c in words[1:]]
        elif words[0]=='coef_culm_theta':
            coef_culm_theta=[float(c) for c in words[1:]]
        elif words[0]=='coef_culm_Rd':
            coef_culm_Rd=[float(c) for c in words[1:]]
        elif words[0]=='coef_culm_Asat_airT':
            coef_culm_Asat_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_culm_Rd_airT':
            coef_culm_Rd_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_grain_Asat':
            coef_grain_Asat=[float(c) for c in words[1:]]
        elif words[0]=='coef_grain_alpha':
            coef_grain_alpha=[float(c) for c in words[1:]]
        elif words[0]=='coef_grain_theta':
            coef_grain_theta=[float(c) for c in words[1:]]
        elif words[0]=='coef_grain_Rd':
            coef_grain_Rd=[float(c) for c in words[1:]]
        elif words[0]=='coef_grain_Asat_airT':
            coef_grain_Asat_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_grain_Rd_airT':
            coef_grain_Rd_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_awn_Asat':
            coef_awn_Asat=[float(c) for c in words[1:]]
        elif words[0]=='coef_awn_alpha':
            coef_awn_alpha=[float(c) for c in words[1:]]
        elif words[0]=='coef_awn_theta':
            coef_awn_theta=[float(c) for c in words[1:]]
        elif words[0]=='coef_awn_Rd':
            coef_awn_Rd=[float(c) for c in words[1:]]
        elif words[0]=='coef_awn_Asat_airT':
            coef_awn_Asat_airT=[float(c) for c in words[1:]]
        elif words[0]=='coef_awn_Rd_airT':
            coef_awn_Rd_airT=[float(c) for c in words[1:]]
            
        if words[0]=='coef_leaf_Asat_coef':
            coef_leaf_Asat_coef=float(words[1])
        elif words[0]=='coef_leaf_alpha_coef':
            coef_leaf_alpha_coef=float(words[1])
        elif words[0]=='coef_leaf_theta_coef':
            coef_leaf_theta_coef=float(words[1])
        elif words[0]=='coef_leaf_Rd_coef':
            coef_leaf_Rd_coef=float(words[1])
        elif words[0]=='coef_culm_Asat_coef':
            coef_culm_Asat_coef=float(words[1])
        elif words[0]=='coef_culm_alpha_coef':
            coef_culm_alpha_coef=float(words[1])
        elif words[0]=='coef_culm_theta_coef':
            coef_culm_theta_coef=float(words[1])
        elif words[0]=='coef_culm_Rd_coef':
            coef_culm_Rd_coef=float(words[1])
        elif words[0]=='coef_grain_Asat_coef':
            coef_grain_Asat_coef=float(words[1])
        elif words[0]=='coef_grain_alpha_coef':
            coef_grain_alpha_coef=float(words[1])
        elif words[0]=='coef_grain_theta_coef':
            coef_grain_theta_coef=float(words[1])
        elif words[0]=='coef_grain_Rd_coef':
            coef_grain_Rd_coef=float(words[1])
        elif words[0]=='coef_awn_Asat_coef':
            coef_awn_Asat_coef=float(words[1])
        elif words[0]=='coef_awn_alpha_coef':
            coef_awn_alpha_coef=float(words[1])
        elif words[0]=='coef_awn_theta_coef':
            coef_awn_theta_coef=float(words[1])
        elif words[0]=='coef_awn_Rd_coef':
            coef_awn_Rd_coef=float(words[1])
            
        elif words[0]=='leaf_NC':
            leaf_NC=[float(c) for c in words[1:]]
        elif words[0]=='grain_NC':
            grain_NC=[float(c) for c in words[1:]]
        elif words[0]=='awn_NC':
            awn_NC=[float(c) for c in words[1:]]
        elif words[0]=='culm_NC':
            culm_NC=[float(c) for c in words[1:]]
    for line in open(config_rayTracing_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='start_time':
            start_time=int(words[1])
        elif words[0]=='end_time':
            end_time=int(words[1])
        elif words[0]=='time_interv':
            time_interv=int(words[1])
        elif words[0]=='lat_field':
            lat_field=float(words[1])
        elif words[0]=='xmin':
            xmin=float(words[1])
        elif words[0]=='xmax':
            xmax=float(words[1])
        elif words[0]=='ymin':
            ymin=float(words[1])
        elif words[0]=='ymax':
            ymax=float(words[1])
        elif words[0]=='zmin':
            zmin=float(words[1])
        elif words[0]=='zmax':
            zmax=float(words[1])
    for line in open(config_canopyPhotosynthesis_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='CP_analysis':
            CP_analysis=int(words[1])
    
    coef_leaf_Asat=np.array(coef_leaf_Asat)*coef_leaf_Asat_coef
    coef_leaf_alpha=np.array(coef_leaf_alpha)*coef_leaf_alpha_coef
    coef_leaf_theta=np.array(coef_leaf_theta)*coef_leaf_theta_coef
    coef_leaf_Rd=np.array(coef_leaf_Rd)*coef_leaf_Rd_coef
    coef_culm_Asat=np.array(coef_culm_Asat)*coef_culm_Asat_coef
    coef_culm_alpha=np.array(coef_culm_alpha)*coef_culm_alpha_coef
    coef_culm_theta=np.array(coef_culm_theta)*coef_culm_theta_coef
    coef_culm_Rd=np.array(coef_culm_Rd)*coef_culm_Rd_coef
    coef_grain_Asat=np.array(coef_grain_Asat)*coef_grain_Asat_coef
    coef_grain_alpha=np.array(coef_grain_alpha)*coef_grain_alpha_coef
    coef_grain_theta=np.array(coef_grain_theta)*coef_grain_theta_coef
    coef_grain_Rd=np.array(coef_grain_Rd)*coef_grain_Rd_coef
    coef_awn_Asat=np.array(coef_awn_Asat)*coef_awn_Asat_coef
    coef_awn_alpha=np.array(coef_awn_alpha)*coef_awn_alpha_coef
    coef_awn_theta=np.array(coef_awn_theta)*coef_awn_theta_coef
    coef_awn_Rd=np.array(coef_awn_Rd)*coef_awn_Rd_coef
    
    DateTimeStr_list,HOUR_list,inDL_list,inSL_list,VPD_list,AT_list,inTL_list=inputFileProcessing(weatherFile)
    given_timePts_num=len(HOUR_list)
    mid_dateTime=datetime.strptime(DateTimeStr_list[given_timePts_num//2], "%Y-%m-%d-%H-%M-%S")
    startDate_of_year=datetime.strptime(DateTimeStr_list[0].split('-')[0]+'-1-1', "%Y-%m-%d")
    current_DOY=(mid_dateTime - startDate_of_year).days + 1
    mid_date_str='-'.join(DateTimeStr_list[given_timePts_num//2].split('-')[0:3])
    fnOutDL=os.path.join(temp_dir,'Run'+str(runID)+'_PPFD_fastTracerDL_'+str(mid_date_str)+'_'+fnOut_remark+'.txt')
    fnOutSL=os.path.join(temp_dir,'Run'+str(runID)+'_PPFD_fastTracerSL_'+str(mid_date_str)+'_'+fnOut_remark+'.txt')
    weatherfnOutDL=weatherFile[0:-4]+'_rayTracingDL.txt'
    weatherfnOutSL=weatherFile[0:-4]+'_rayTracingSL.txt'
    weatherOut_DL=open(weatherfnOutDL,'w')
    weatherOut_SL=open(weatherfnOutSL,'w')
    YEAR=DateTimeStr_list[0].split('-')[0]
    for i in range(start_time,end_time+1,time_interv):
        content=[YEAR,current_DOY,i-0.01,25,80,1000,0] # year, DOY, time, Temp, RH, TotPPFD, ScatterPPFD
        weatherOut_DL.write('\t'.join([str(c) for c in content])+'\n')
        content=[YEAR,current_DOY,i,25,80,1000,0] # year, DOY, time, Temp, RH, TotPPFD, ScatterPPFD
        weatherOut_DL.write('\t'.join([str(c) for c in content])+'\n')
        content=[YEAR,current_DOY,i+0.01,25,80,1000,0] # year, DOY, time, Temp, RH, TotPPFD, ScatterPPFD
        weatherOut_DL.write('\t'.join([str(c) for c in content])+'\n')
        
    content=[YEAR,current_DOY,6,25,80,1000,1000] # year, DOY, time, Temp, RH, TotPPFD, ScatterPPFD
    weatherOut_SL.write('\t'.join([str(c) for c in content])+'\n')
    weatherOut_DL.close()
    weatherOut_SL.close()
    
    coef_leaf_AQ_paras=[coef_leaf_Asat,coef_leaf_alpha,coef_leaf_theta,coef_leaf_Rd,coef_leaf_Asat_airT,coef_leaf_Rd_airT]
    coef_culm_AQ_paras=[coef_culm_Asat,coef_culm_alpha,coef_culm_theta,coef_culm_Rd,coef_culm_Asat_airT,coef_culm_Rd_airT]
    coef_grain_AQ_paras=[coef_grain_Asat,coef_grain_alpha,coef_grain_theta,coef_grain_Rd,coef_grain_Asat_airT,coef_grain_Rd_airT]
    coef_awn_AQ_paras=[coef_awn_Asat,coef_awn_alpha,coef_awn_theta,coef_awn_Rd,coef_awn_Asat_airT,coef_awn_Rd_airT]
    
    DL_timePts=int((end_time-start_time)/time_interv+1)
    format_DLfn='%d\t'*5+'%.2f\t'*9+'%.2f\t%.3f\t%.3f\t%.4f\t'+'%d\t'*DL_timePts*7
    format_SLfn='%d\t'*5+'%.2f\t'*9+'%.2f\t%.3f\t%.3f\t%.4f\t'+'%d\t'*7
    
    canopy_rayTracing_fn=os.path.join(resultDir,'Run'+str(runID)+'_canopyRay_'+str(mid_date_str)+'_'+fnOut_remark+'.txt')
    content1='plantID\ttillerID\tposID\torganCode\textraID\tx1\ty1\tz1\tx2\ty2\tz2\tx3\ty3\tz3\t'
    content2='organNC\tKt\tKr\torganArea\t'
    content3=''
    for time_temp in HOUR_list:
        content_more=(str(time_temp)+'\t').join(['adDL','adSL','adDFL','abDL','abSL','abDFL','TL',''])
        content3=content3+content_more
    header_str_canopyRay=content1+content2+content3[0:-1]+'\n'
    format_canopyRay='%d\t'*5+'%.2f\t'*9+'%.2f\t%.3f\t%.3f\t%.4f\t'+'%d\t'*given_timePts_num*7
    format_canopyRay=format_canopyRay[0:-1]+'\n'
    
    canopy_photosynthesis_fn=os.path.join(resultDir,'Run'+str(runID)+'_canopyPhotoRate_'+str(mid_date_str)+'_'+fnOut_remark+'.txt')
    content3=''
    for time_temp in HOUR_list:
        content3=content3+str(time_temp)+'Anet\t'
    header_str_canopyPhotoRate=header_str_canopyRay[0:-1]+'\t'+content3[0:-1]+'\n'
    format_canopyPhotoRate=format_canopyRay[0:-1]+'\t'+'%.2f\t'*given_timePts_num # umol/m2/s
    format_canopyPhotoRate=format_canopyPhotoRate[0:-1]+'\n'
    
    if not os.path.exists(canopy_photosynthesis_fn):
        if not os.path.exists(canopy_rayTracing_fn):
            if not (os.path.exists(fnOutDL) and os.path.exists(fnOutSL)):
                if not os.path.exists(M_plant_clean_fn):
                    print('Canopy structure not exists! Canopy structure reconstructing ...\n')
                    canopyReconstructor(canopyConstructorLoc,inputDir,keyPoints_raw_fn_short,M_plant_fn,M_plant_clean_fn,img_H,cultivarID,runID,config_morphology_fn,config_physiology_fn,config_rayTracing_fn,overWrite,withSpike)
                #print('Canopy ray tracing ...\n')
                fieldRayTracer(M_plant_clean_fn,config_rayTracing_fn,weatherfnOutDL,weatherfnOutSL,fnOutDL,fnOutSL,fastTracerLoc,current_DOY,DL_timePts,format_DLfn,format_SLfn)
            else:
                print('Existing standard DL/SL ray tracing results found!\n')
            DL_onOrgan_time,SL_onOrgan_time,canopy_mat,facetNum=DLSL_extractor(fnOutDL,fnOutSL,weatherFile,config_rayTracing_fn)
            totLight=DL_onOrgan_time+SL_onOrgan_time # totLight is absorbed light but not incident light
            canopy_rayTracing_mat=np.hstack((canopy_mat,totLight))
            # plot figures and generate video
            #rayTracing_dir=os.path.join(resultDir,'rayTracingFigs')
            #if not os.path.exists(rayTracing_dir):
            #    os.mkdir(rayTracing_dir)
            #print_info=print('Deleting temporary and outdated files ...');
            #delete(fnOutDL)
            #delete(fnOutSL)
            # save to file
            print('Canopy ray tracing done. Writing canopy ray tracing result into file ...\n')
            saveData(canopy_rayTracing_fn,header_str_canopyRay,canopy_rayTracing_mat,format_canopyRay,facetNum)
            #if os.path.exists(canopy_stucture_mat_fn):
            #    delete(canopy_stucture_mat_fn)
        else:
            print('Existing canopy ray tracing results found! Loading it ...\n')
            facetNum=fileLineNumFastCounter(canopy_rayTracing_fn)-1
            colNum=18+given_timePts_num*7
            canopy_rayTracing_mat=np.zeros([facetNum,colNum])
            line_i=0
            with open(canopy_rayTracing_fn,'r') as f:
                f.readline()
                for line in f:
                    canopy_rayTracing_mat[line_i]=[float(c) for c in line.strip().split('\t')]
                    line_i+=1
        print('Simulating canopy photosynthesis ...\n')  
        organPos =canopy_rayTracing_mat[:,2] # soil:-6, awn:-5, glume:-4, grain:-3, branch:-2, sheath:-1, stem:0, leaf:1..10 
        #organCode=canopy_rayTracing_mat[:,3] # soil:1, stem:3, leaf:4, spike:5 
        organNC=canopy_rayTracing_mat[:,14]
        organKt=canopy_rayTracing_mat[:,15]
        organKr=canopy_rayTracing_mat[:,16]
        facet_area=canopy_rayTracing_mat[:,17]
        totLight=canopy_rayTracing_mat[:,17+7:17+7*(given_timePts_num+1):7] # only extract TotPPFD
        
        A_mat=organPhotoCalculator(organPos,facet_area, organKt,organKr,organNC, totLight, coef_leaf_AQ_paras,coef_culm_AQ_paras,coef_grain_AQ_paras,coef_awn_AQ_paras,facetNum,AT_list) # umol/m2/s
        
        canopy_photosynthesis_mat=np.hstack((canopy_rayTracing_mat,A_mat)) # umol/m2/s
        # save to file
        print('Canopy photosynthesis calculation done. Writing canopy photosynthesis result into file ...\n')
        saveData(canopy_photosynthesis_fn,header_str_canopyPhotoRate,canopy_photosynthesis_mat,format_canopyPhotoRate,facetNum)
        print('Deleting temporary and outdated files ...')
        try:
            os.remove(canopy_rayTracing_fn)
        except:
            pass
        try:
            os.remove(M_plant_clean_fn)
        except:
            pass
    else:
        print('Existing canopy photosynthesis result found! Loading it ...\n')
        facetNum=fileLineNumFastCounter(canopy_photosynthesis_fn)-1
        colNum=18+given_timePts_num*(7+1)
        canopy_photosynthesis_mat=np.zeros([facetNum,colNum])
        line_i=0
        with open(canopy_photosynthesis_fn,'r') as f:
            f.readline()
            for line in f:
                canopy_photosynthesis_mat[line_i]=[float(c) for c in line.strip().split('\t')]
                line_i+=1
    # canopy photosynthesis analysis
    if CP_analysis:
        print('Plotting and saving canopy photosynthesis analytical results ...\n')
        canopy_photosynthesis_mat=canopy_photosynthesis_mat[canopy_photosynthesis_mat[:,3]>1,:] # exclude soil
        plot_canopyPhotoStats(canopy_photosynthesis_mat,DateTimeStr_list,HOUR_list,inDL_list,inSL_list,VPD_list,AT_list,inTL_list,xmin,xmax,ymin,ymax,resultDir,mid_date_str,runID,fnOut_remark)
    
    eclipse=time.time()-start
    print('All done! Time used: %.2fs.\n'%eclipse)
    
    