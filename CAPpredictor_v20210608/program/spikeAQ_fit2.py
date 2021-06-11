#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
spikeAQ_fit2.py

1. construct 3D wheat spike with and without awns
2. perform ray tracing simulation with and without awns
3. fit AQ curve with and without awns
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

import matplotlib.pyplot as plt
#from scipy.interpolate import pchip
from scipy.interpolate import PchipInterpolator
from scipy.optimize import leastsq,curve_fit,fmin_slsqp


global epsilon
epsilon=0.0001
global rayTrace_fit_error_min
rayTrace_fit_error_min=0.1

# interpolant = PchipInterpolator(x, y)
# xnew = np.linspace(0, 10, num=41, endpoint=True)
# interpolant(xnew)

__author__="CHANG Tiangen"
__date__="20201207"
__copyright__="Copyright 2020 Aspar CHANG"
__license__="SIPPE"
__version__="0.1"
__email__="changtiangen@sippe.ac.cn"
__status__="Prototype"
__credits__=[]

    
def createHelp():
    """
    Create the command line interface of the programme.
    """
  
    epilog_string="Any bug is welcome reported to changtiangen@sippe.ac.cn"
    description_string='The program is going to fit the wheat AQ curve.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-3d', '--3dfnOut', dest='file_3d', help='3d output file')
    parser.add_argument('-i', '--fnIn2', dest='file_in2', help='Input spike basic info file')
    parser.add_argument('-iw', '--fnIn_w', dest='file_in_w', help='Input weather data file')
    parser.add_argument('-o', '--fnOut', dest='file_out', help='Output file')
    parser.add_argument('-ih', '--imageH', dest='img_H', type=int, default=3264, help='Image height in pixels')
    parser.add_argument('-lineNA', '--cultivar', dest='lineNA', type=str, default='', help='Cultivar specific key word in image file name')
    parser.add_argument('-spikeAng', '--spikeAng', dest='spikeAng', type=int,nargs='*', default=[0], help='Ear rotation angles in P chamber')
    parser.add_argument('-rID', '--runID', dest='runID', type=int, default=1, help='Remark of current run ID, used as file prefix and random seed.')
    parser.add_argument('-ow', '--overWrite', dest='ow', type=int, default=1, help='Overwrite the existing files (1) or not (0).')
    parser.add_argument('-rd', '--rayDense', dest='RD', type=float, default=0.05, help='Ray density (cm).')
    parser.add_argument('-show', '--show', dest='showOn', type=int, default=0, help='Show 3D structure on screen (1) or not (0).')
    op=parser.parse_args()
    return op

def rayFileMerge(M_plant_ray_fn_list,M_plant_ray_fn,ave_num):
    data1_raw=open(M_plant_ray_fn_list[0],'r').readlines()
    col_num=len(data1_raw[0].strip().split('\t'))
    row_num=len(data1_raw)
    data_clean=np.zeros([row_num-1,col_num])
    for i in range(1,row_num):
        line=data1_raw[i].replace('-nan(ind)','0')
        words=line.strip().split('\t')
        if len(words)>1:
            data_clean[i-1]=[float(c) for c in words]
    for file_i in range(1,len(M_plant_ray_fn_list)):
        data1_raw=open(M_plant_ray_fn_list[file_i],'r').readlines()
        for i in range(1,row_num):
            line=data1_raw[i].replace('-nan(ind)','0')
            words=line.strip().split('\t')
            if len(words)>1:
                PPFDs=[float(c) for c in words[18:]]
                data_clean[i-1][18:]=data_clean[i-1][18:]+PPFDs
    if ave_num>1:
        data_clean[:,18:]=data_clean[:,18:]/ave_num
    fnOut=open(M_plant_ray_fn,'w')
    fnOut.write(data1_raw[0])
    for i in range(row_num-1):
        content1='\t'.join([str(round(c)) for c in data_clean[i][0:5]])
        content2='\t'.join([str(c) for c in data_clean[i][5:]])
        content=content1+'\t'+content2+'\n'
        fnOut.write(content)
    fnOut.close()
    return data_clean

def uniqueFitJudger(itera_test_in,fit_error_new,fit_error,popt,popt_old,uniq_old):
    itera_test_out=itera_test_in
    fit_error_out=fit_error
    popt_new=popt_old
    theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=popt_old
    theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq=uniq_old
    if fit_error_new<fit_error-epsilon:
        print('fit_error_new: %.2f'%fit_error_new)
        popt_new=popt
        fit_error_out=fit_error_new
        itera_test_out=1
        theta_uniq='YES'
        Amax25_uniq='YES'
        phi_CO2_uniq='YES'
        Rd25_uniq='YES'
        theta_awn_uniq='YES'
        Amax25_awn_uniq='YES'
        phi_CO2_awn_uniq='YES'
        Rd25_awn_uniq='YES'
    elif fit_error_new>fit_error+epsilon:
        True
    else:
        if abs(theta-popt[0])/(theta+epsilon)>0.05:
            theta_uniq='NO'
        if abs(Amax25-popt[1])/(Amax25+epsilon)>0.05:
            Amax25_uniq='NO'
        if abs(phi_CO2-popt[2])/(phi_CO2+epsilon)>0.05:
            phi_CO2_uniq='NO'
        if abs(Rd25-popt[3])/(Rd25+epsilon)>0.05:
            Rd25_uniq='NO'
        if abs(theta_awn-popt[4])/(theta_awn+epsilon)>0.05:
            theta_awn_uniq='NO'
        if abs(Amax25_awn-popt[5])/(Amax25_awn+epsilon)>0.05:
            Amax25_awn_uniq='NO'
        if abs(phi_CO2_awn-popt[6])/(phi_CO2_awn+epsilon)>0.05:
            phi_CO2_awn_uniq='NO'
        if abs(Rd25_awn-popt[7])/(Rd25_awn+epsilon)>0.05:
            Rd25_awn_uniq='NO'
    uniq_new=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq]
    return itera_test_out,fit_error_out,popt_new,uniq_new

def fit_error_reset(fit_error):
    if math.isnan(fit_error):
        fit_error_out=100000000
    else:
        fit_error_out=fit_error
    return fit_error_out

def fit_AQ(grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,PAR_ave,A_ave,CT_ave,A_noAwn_ave,CT_noAwn_ave,Amax_oef,Rd_coef):
    
    x=np.array(PAR_ave)
    y=np.array(A_ave)
    z=np.array(CT_ave)
    y2=np.array(A_noAwn_ave)
    z2=np.array(CT_noAwn_ave)
    
    tot_PPF_grain1=np.dot(grain1_PPFD_array,grain1_area_array)/10 # nmol s-1
    x_tot_grain1=np.array(PAR_ave)/2000*tot_PPF_grain1
    tot_PPF_awn1=np.dot(awn1_PPFD_array,awn1_area_array)/10 # nmol s-1
    x_tot_awn1=np.array(PAR_ave)/2000*tot_PPF_awn1
    tot_PPF_grain2=np.dot(grain2_PPFD_array,grain2_area_array)/10 # nmol s-1
    x_tot_grain2=np.array(PAR_ave)/2000*tot_PPF_grain2
    
    ################ simple total absorbed PPFD to A AQ fitting ################
    itera_test=0
    fit_error=100000000
    theta=0
    Amax25=0
    phi_CO2=0
    Rd25=0
    theta_uniq='YES'
    Amax25_uniq='YES'
    phi_CO2_uniq='YES'
    Rd25_uniq='YES'
    theta_awn=0
    Amax25_awn=0
    phi_CO2_awn=0
    Rd25_awn=0
    theta_awn_uniq='YES'
    Amax25_awn_uniq='YES'
    phi_CO2_awn_uniq='YES'
    Rd25_awn_uniq='YES'
    
    while itera_test < 50:
        print('            Simple AQ fitting: iteration %d, fit_error %.1f'%(itera_test,fit_error))
        theta_ini=random.uniform(0.01,0.999)
        Amax25_ini=random.randint(1,100)*1.01
        phi_CO2_ini=random.uniform(0,0.1)
        Rd25_ini=random.randint(1,100)*1.01
        theta_awn_ini=random.uniform(0.01,0.999)
        Amax25_awn_ini=random.randint(1,100)*1.01
        phi_CO2_awn_ini=random.uniform(0,0.1)
        Rd25_awn_ini=random.randint(1,100)*1.01
        paras0=np.array([theta_ini,Amax25_ini,phi_CO2_ini,Rd25_ini,theta_awn_ini,Amax25_awn_ini,phi_CO2_awn_ini,Rd25_awn_ini]) # theta,Amax25,alpha,Rd
        popt= fmin_slsqp(fit_residuals, paras0, bounds=[(0.01,0.999),(0,101),(0,0.1),(0,101),(0.01,0.999),(0,101),(0,0.1),(0,101)], args=(x_tot_grain1,x_tot_awn1,y,z,x_tot_grain2,y2,z2,Amax_coef,Rd_coef),iprint = 0)
        popt=popt[0:]
        if math.isnan(popt[0]):
            continue
        if theta+Amax25+phi_CO2+theta_awn+Amax25_awn+phi_CO2_awn==0:
            theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=popt
            fit_error=fit_residuals(popt[0:],x_tot_grain1,x_tot_awn1,y,z,x_tot_grain2,y2,z2,Amax_coef,Rd_coef)
            #fit_error=fit_error_reset(fit_error)
        else:
            fit_error_new=fit_residuals(popt[0:],x_tot_grain1,x_tot_awn1,y,z,x_tot_grain2,y2,z2,Amax_coef,Rd_coef)
            #fit_error_new=fit_error_reset(fit_error_new)
            popt_old=[theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn]
            uniq_old=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq]
            itera_test,fit_error,popt_new,uniq_new = uniqueFitJudger(itera_test,fit_error_new,fit_error,popt,popt_old,uniq_old)
            theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=popt_new
            theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq=uniq_new
        itera_test=itera_test+1
    Rd_array=Rd25*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Rd_awn_array=Rd25_awn*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Amax_array=Amax25*(Amax_coef[0]*z**2 + Amax_coef[1]*z + Amax_coef[2]) # Rd-T empirical relation
    Amax_awn_array=Amax25_awn*(Amax_coef[0]*z**2 + Amax_coef[1]*z + Amax_coef[2]) # Rd-T empirical relation
    Ag_max=Amax_array+Rd_array
    Ag_awn_max=Amax_awn_array+Rd_awn_array
    A_grain1_fit_array=(phi_CO2*x_tot_grain1+Ag_max-np.sqrt((phi_CO2*x_tot_grain1+Ag_max)**2-4*theta*phi_CO2*x_tot_grain1*Ag_max))/(2*theta)-Rd_array
    A_awn1_fit_array=(phi_CO2_awn*x_tot_awn1+Ag_awn_max-np.sqrt((phi_CO2_awn*x_tot_awn1+Ag_awn_max)**2-4*theta_awn*phi_CO2_awn*x_tot_awn1*Ag_awn_max))/(2*theta_awn)-Rd_awn_array
    Rd_grain2_array=Rd25*(Rd_coef[0]*z2**2 + Rd_coef[1]*z2 + Rd_coef[2]) # Rd-T empirical relation
    Amax_grain2_array=Amax25*(Amax_coef[0]*z2**2 + Amax_coef[1]*z2 + Amax_coef[2]) # Rd-T empirical relation
    
    Ag_max=Amax_grain2_array+Rd_grain2_array
    A_fit_grain2_array=(phi_CO2*x_tot_grain2+Ag_max-np.sqrt((phi_CO2*x_tot_grain2+Ag_max)**2-4*theta*phi_CO2*x_tot_grain2*Ag_max))/(2*theta)-Rd_grain2_array
    
    A_fit_all = np.hstack((A_grain1_fit_array + A_awn1_fit_array,A_fit_grain2_array))
    y_entire=np.hstack((y,y2))
    R2_fit=1-sum((A_fit_all-y_entire)**2)/sum((y_entire-np.mean(y_entire))**2)
    result_1=[theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq,fit_error,R2_fit,A_fit_all]
    
    ################ Ray tracing PAR_on_facet VS A_ave AQ fitting ################
    itera_test=0
    fit_error=100000000
    theta=0
    Amax25=0
    phi_CO2=0
    Rd25=0
    theta_uniq='YES'
    Amax25_uniq='YES'
    phi_CO2_uniq='YES'
    Rd25_uniq='YES'
    theta_awn=0
    Amax25_awn=0
    phi_CO2_awn=0
    Rd25_awn=0
    theta_awn_uniq='YES'
    Amax25_awn_uniq='YES'
    phi_CO2_awn_uniq='YES'
    Rd25_awn_uniq='YES'
    Amax25_ub=20
    Rd25_ub=10
    Amax25_awn_ub=20
    Rd25_awn_ub=10
    while itera_test < 50:
        print('            Ray tracing AQ fitting: iteration %d, fit_error %.1f'%(itera_test,fit_error))
        theta_ini=random.uniform(0.01,0.999)
        Amax25_ini=random.uniform(0,Amax25_ub)
        phi_CO2_ini=random.uniform(0,0.1)
        Rd25_ini=random.uniform(0,Rd25_ub)
        theta_awn_ini=random.uniform(0.01,0.999)
        Amax25_awn_ini=random.uniform(0,Amax25_awn_ub)
        phi_CO2_awn_ini=random.uniform(0,0.1)
        Rd25_awn_ini=random.uniform(0,Rd25_awn_ub)
#         theta_ini=0.839100165
#         Amax25_ini=4.093976743
#         phi_CO2_ini=0.079715078
#         Rd25_ini=2.420171864
#         theta_awn_ini=1
#         Amax25_awn_ini=9.406281688
#         phi_CO2_awn_ini=0.1
#         Rd25_awn_ini=0.1
        
        paras0=np.array([theta_ini,Amax25_ini,phi_CO2_ini,Rd25_ini,theta_awn_ini,Amax25_awn_ini,phi_CO2_awn_ini,Rd25_awn_ini]) # theta,Amax25,alpha,Rd
        popt= fmin_slsqp(fit_residuals2, paras0, bounds=[(0.01,0.999),(0,Amax25_ub),(0,0.1),(0,Rd25_ub),(0.01,0.999),(0,Amax25_awn_ub),(0,0.1),(0,Rd25_awn_ub)], args=(grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,x,y,z,y2,z2,Amax_coef,Rd_coef),iprint = 0)
        popt=popt[0:]
        if math.isnan(popt[0]):
            continue
        if theta+Amax25+phi_CO2+theta_awn+Amax25_awn+phi_CO2_awn==0:
            theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=popt
            fit_error=fit_residuals2(popt[0:],grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,x,y,z,y2,z2,Amax_coef,Rd_coef)
        else:
            fit_error_new=fit_residuals2(popt[0:],grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,x,y,z,y2,z2,Amax_coef,Rd_coef)
            popt_old=[theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn]
            uniq_old=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq]
            itera_test,fit_error,popt_new,uniq_new = uniqueFitJudger(itera_test,fit_error_new,fit_error,popt,popt_old,uniq_old)
            theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=popt_new
            theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq=uniq_new
        if fit_error<rayTrace_fit_error_min:
            break
        itera_test=itera_test+1
    AQ_paras=[theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn]
    A_fit_array2=Aspike_fit(grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,AQ_paras,x,z,z2,Amax_coef,Rd_coef)
    R2_fit=1-sum((A_fit_array2-y_entire)**2)/sum((y_entire-np.mean(y_entire))**2)
    result_2=[theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,theta_awn_uniq,Amax25_awn_uniq,phi_CO2_awn_uniq,Rd25_awn_uniq,fit_error,R2_fit,A_fit_array2]

    return x_tot_grain1,x_tot_awn1,x_tot_grain2,result_1,result_2
    
def Aspike_fit(grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,AQ_paras,PAR_ave,CT_ave,CT_noAwn_ave,Amax_coef,Rd_coef):
    theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=AQ_paras
    y_fit=[]
    for PAR_i in range(len(PAR_ave)):
        x=grain1_PPFD_array*PAR_ave[PAR_i]/2000
        Rd=Rd25*(Rd_coef[0]*CT_ave[PAR_i]**2 + Rd_coef[1]*CT_ave[PAR_i] + Rd_coef[2]) # Rd-T empirical relation
        Amax=Amax25*(Amax_coef[0]*CT_ave[PAR_i]**2 + Amax_coef[1]*CT_ave[PAR_i] + Amax_coef[2]) # Rd-T empirical relation
        Ag_max=Amax+Rd
        A_fit_grain1_array=(phi_CO2*x+Ag_max-np.sqrt((phi_CO2*x+Ag_max)**2-4*theta*phi_CO2*x*Ag_max))/(2*theta)-Rd
        x=awn1_PPFD_array*PAR_ave[PAR_i]/2000
        Rd_awn1=Rd25_awn*(Rd_coef[0]*CT_ave[PAR_i]**2 + Rd_coef[1]*CT_ave[PAR_i] + Rd_coef[2]) # Rd-T empirical relation
        Amax_awn1=Amax25_awn*(Amax_coef[0]*CT_ave[PAR_i]**2 + Amax_coef[1]*CT_ave[PAR_i] + Amax_coef[2]) # Rd-T empirical relation

        Ag_max_awn=Amax_awn1+Rd_awn1
        A_fit_awn1_array=(phi_CO2_awn*x+Ag_max_awn-np.sqrt((phi_CO2_awn*x+Ag_max_awn)**2-4*theta_awn*phi_CO2_awn*x*Ag_max_awn))/(2*theta_awn)-Rd_awn1
        y_fit.append(float(np.dot(grain1_area_array,A_fit_grain1_array)/10+np.dot(awn1_area_array,A_fit_awn1_array)/10)) # nmol ear-1 s-1
    for PAR_i in range(len(PAR_ave)):
        x=grain2_PPFD_array*PAR_ave[PAR_i]/2000
        Rd=Rd25*(Rd_coef[0]*CT_noAwn_ave[PAR_i]**2 + Rd_coef[1]*CT_noAwn_ave[PAR_i] + Rd_coef[2]) # Rd-T empirical relation
        Amax=Amax25*(Amax_coef[0]*CT_noAwn_ave[PAR_i]**2 + Amax_coef[1]*CT_noAwn_ave[PAR_i] + Amax_coef[2]) # Rd-T empirical relation
        Ag_max=Amax+Rd
        A_fit_grain2_array=(phi_CO2*x+Ag_max-np.sqrt((phi_CO2*x+Ag_max)**2-4*theta*phi_CO2*x*Ag_max))/(2*theta)-Rd
        y_fit.append(float(np.dot(grain2_area_array,A_fit_grain2_array)/10)) # nmol ear-1 s-1
    return np.array(y_fit)
    
def fit_residuals(paras,x_grain1,x_awn1,y,z,x_grain2,y2,z2,Amax_coef,Rd_coef):#y=A_entire, x=PPFD, z=CT, y2=A_noAwn, z2=CT2
    theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=paras
    Rd_grain1_array=Rd25*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Rd_awn1_array=Rd25_awn*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Rd_grain2_array=Rd25*(Rd_coef[0]*z2**2 + Rd_coef[1]*z2 + Rd_coef[2]) # Rd-T empirical relation
    Amax_grain1_array=Amax25*(Amax_coef[0]*z**2 + Amax_coef[1]*z + Amax_coef[2]) # Amax-T empirical relation
    Amax_awn1_array=Amax25_awn*(Amax_coef[0]*z**2 + Amax_coef[1]*z + Amax_coef[2]) # Amax-T empirical relation
    Amax_grain2_array=Amax25*(Amax_coef[0]*z2**2 + Amax_coef[1]*z2 + Amax_coef[2]) # Amax-T empirical relation
    Ag_grain1_max=Amax_grain1_array+Rd_grain1_array
    Ag_awn1_max=Amax_awn1_array+Rd_awn1_array
    Ag_grain2_max=Amax_grain2_array+Rd_grain2_array
    
    sqrt1=(phi_CO2*x_grain1+Ag_grain1_max)**2-4*theta*phi_CO2*x_grain1*Ag_grain1_max
    sqrt2=(phi_CO2_awn*x_awn1+Ag_awn1_max)**2-4*theta_awn*phi_CO2_awn*x_awn1*Ag_awn1_max
    sqrt3=(phi_CO2*x_grain2+Ag_grain2_max)**2-4*theta*phi_CO2*x_grain2*Ag_grain2_max
    y_grain1_fit=(phi_CO2*x_grain1+Ag_grain1_max-np.sqrt(sqrt1))/(2*theta)-Rd_grain1_array
    y_awn1_fit=(phi_CO2_awn*x_awn1+Ag_awn1_max-np.sqrt(sqrt2))/(2*theta_awn)-Rd_awn1_array
    y_fit1=y_grain1_fit+y_awn1_fit
    y_grain2_fit=(phi_CO2*x_grain2+Ag_grain2_max-np.sqrt(sqrt3))/(2*theta)-Rd_grain2_array
    y_fit=np.hstack((y_fit1,y_grain2_fit))
    y_entire=np.hstack((y,y2))
    fit_err=sum((y_fit-y_entire)**2)
    #print('Simple fit: ',y_fit,y_entire)
    return fit_err
    
def fit_residuals2(paras,grain1_area_array,grain1_PPFD_array,awn1_area_array,awn1_PPFD_array,grain2_area_array,grain2_PPFD_array,PAR_ave,A_ave,CT_ave,A_noAwn_ave,CT_noAwn_ave,Amax_coef,Rd_coef):
    theta,Amax25,phi_CO2,Rd25,theta_awn,Amax25_awn,phi_CO2_awn,Rd25_awn=paras
    Rd_grain1_array=Rd25*(Rd_coef[0]*CT_ave**2 + Rd_coef[1]*CT_ave + Rd_coef[2]) # Rd-T empirical relation
    Rd_awn1_array=Rd25_awn*(Rd_coef[0]*CT_ave**2 + Rd_coef[1]*CT_ave + Rd_coef[2]) # Rd-T empirical relation
    Rd_grain2_array=Rd25*(Rd_coef[0]*CT_noAwn_ave**2 + Rd_coef[1]*CT_noAwn_ave + Rd_coef[2]) # Rd-T empirical relation
    Amax_grain1_array=Amax25*(Amax_coef[0]*CT_ave**2 + Amax_coef[1]*CT_ave + Amax_coef[2]) # Amax-T empirical relation
    Amax_awn1_array=Amax25_awn*(Amax_coef[0]*CT_ave**2 + Amax_coef[1]*CT_ave + Amax_coef[2]) # Amax-T empirical relation
    Amax_grain2_array=Amax25*(Amax_coef[0]*CT_noAwn_ave**2 + Amax_coef[1]*CT_noAwn_ave + Amax_coef[2]) # Amax-T empirical relation
    
    y_fit=[]
    for PAR_i in range(len(PAR_ave)):
        x_grain1=grain1_PPFD_array*(PAR_ave[PAR_i]/2000)
        x_awn1=awn1_PPFD_array*(PAR_ave[PAR_i]/2000)
        Rd=Rd_grain1_array[PAR_i]
        Rd_awn1=Rd_awn1_array[PAR_i]
        Amax=Amax_grain1_array[PAR_i]
        Amax_awn1=Amax_awn1_array[PAR_i]
        Ag_max=Amax+Rd
        Ag_max_awn=Amax_awn1+Rd_awn1
        
        sqrt_grain1=(phi_CO2*x_grain1+Ag_max)**2-(4*theta*phi_CO2*Ag_max)*x_grain1
        sqrt_awn1=(phi_CO2_awn*x_awn1+Ag_max_awn)**2-(4*theta_awn*phi_CO2_awn*Ag_max_awn)*x_awn1
        A_fit_grain1_array=(phi_CO2*x_grain1+Ag_max-np.sqrt(sqrt_grain1))/(2*theta)-Rd
        A_fit_awn1_array=(phi_CO2_awn*x_awn1+Ag_max_awn-np.sqrt(sqrt_awn1))/(2*theta_awn)-Rd_awn1
        y_fit.append(float(np.dot(grain1_area_array,A_fit_grain1_array)/10+np.dot(awn1_area_array,A_fit_awn1_array)/10)) # nmol ear-1 s-1
    for PAR_i in range(len(PAR_ave)):
        x=grain2_PPFD_array*(PAR_ave[PAR_i]/2000)
        Rd=Rd_grain2_array[PAR_i]
        Amax=Amax_grain2_array[PAR_i]
        Ag_max=Amax+Rd
        
        sqrt_grain2=(phi_CO2*x+Ag_max)**2-(4*theta*phi_CO2*Ag_max)*x
        A_fit_grain2_array=(phi_CO2*x+Ag_max-np.sqrt(sqrt_grain2))/(2*theta)-Rd
        y_fit.append(float(np.dot(grain2_area_array,A_fit_grain2_array)/10)) # nmol ear-1 s-1
    A_entire=np.hstack((A_ave,A_noAwn_ave))
    #print('rayTrace fit: ',y_fit,A_entire)
    fit_err=sum((np.array(y_fit)-A_entire)**2)
    return fit_err
    
def canopyReconstructor(M_plant_raw_fn,img_H,cultivarID,runID,spikeAng,awn,fnIn_spikeInfo,overWrite):
    # Canopy structure reconstructing
    print('    Spike structure reconstructing ...\n')
    if cultivarID:
        command='python .\wheatSpikeConstructor.py -i %s -o %s -ih %d -lineNA %s -rID %d -ow %d -ang %d -awn %d -log 0'%(fnIn_spikeInfo,M_plant_raw_fn,img_H,cultivarID,runID,overWrite,spikeAng,awn)
    else:
        command='python .\wheatSpikeConstructor.py -i %s -o %s -ih %d -rID %d -ow %d -ang %d -awn %d -log 0'%(fnIn_spikeInfo,M_plant_raw_fn,img_H,runID,overWrite,spikeAng,awn)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)
    
def plantDrawer(M_plant_fn,plant3D_fn,showOn):
    print('    Spike structure plotting ...\n')
    command='python .\plantDraw3D.py -i %s -o %s -ccs 4 -cql 0 -cqu 1 -cm plasma -alpha 1 -va1 60 -show %d'%(M_plant_fn,plant3D_fn,showOn)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command) 
    
def PChamberRayTracer(rayDense,M_plant2_fn,M_plant_ray2_fn,M_plant1_fn,M_plant_ray1_fn,fn_weather,DOY):
    print('P-Chamber spike ray tracing. Face A ...\n')
    command='.\\fastTracerVS2019_0.exe -L 31 -d %d -D -0.01 29.01 -2.01 2.01 0 5.02 -S 12 -n %f -W 6 2 7 -m %s -o %s -C %s -s 1'%(DOY,rayDense,M_plant1_fn,M_plant_ray1_fn,fn_weather)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)   
    print('P-Chamber spike ray tracing. Face B ...\n')
    command='.\\fastTracerVS2019_0.exe -L 31 -d %d -D -0.01 29.01 -2.01 2.01 0 5.02 -S 12 -n %f -W 6 2 7 -m %s -o %s -C %s -s 1'%(DOY,rayDense,M_plant2_fn,M_plant_ray2_fn,fn_weather)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)   
    
def fnIn_fnOut_preprocessing(fnOut_str,fnIn_str,runID):
    fnOut_list=fnOut_str.split('\\')
    resultDir='\\'.join(fnOut_list[0:-1])
    if not resultDir:
        resultDir='.'
    if not os.path.exists(resultDir):
        print('  Result directory not exists, creating it ...')
        try:
            os.mkdir(resultDir)
        except:
            raise Exception('Cannot create result directory: %s'%resultDir)
    if 'Run' not in fnOut_list[-1]:
        fnOut_list[-1]='Run'+str(runID)+'_'+fnOut_list[-1]
    fnOut='\\'.join(fnOut_list)
    fnIn_list=fnIn_str.split('\\')
    if 'Run' not in fnIn_list[-1]:
        fnIn_list[-1]='Run'+str(runID)+'_'+fnIn_list[-1]
    fnIn='\\'.join(fnIn_list)
    return fnOut,fnIn

if __name__=="__main__":
    start=time.time()
    op=createHelp()
    width_format=94
    welcome_info='Welcome to use spike AQ fitter.'
    time_info=time.asctime( time.localtime(time.time()) )
    print(welcome_info.center(width_format, '-'))
    print(time_info.center(width_format, ' '))
    
    runID=op.runID
    random.seed(runID)
    fn_weather=op.file_in_w
    fnIn_spikeInfo=op.file_in2
    img_H=op.img_H
    cultivarID=op.lineNA
    fnOut_str=op.file_out
    overWrite=op.ow
    rayDense=op.RD
    showOn=op.showOn
    rotateAng_list=op.spikeAng # degree
    
    fnOut,op.file_3d=fnIn_fnOut_preprocessing(fnOut_str,op.file_3d,runID)
    op.file_out=fnOut
    data = json.dumps(vars(op), indent=4,ensure_ascii=False, sort_keys=False,separators=(',', ':'))
    print('Running parameters:')
    print(data)
    
    training_data='morning' # 'morning','afternoon'
    for line in open(fnIn_spikeInfo,'r').readlines():
        words=line.strip().split('\t')
        if words[0]=='PPFD':
            PPFD_list=[int(c) for c in words[1:]]
        elif words[0]=='coef_Rd_airT':
            Rd_coef=[float(c) for c in words[1:]] 
        elif words[0]=='coef_Asat_airT':
            Amax_coef=[float(c) for c in words[1:]] 
        elif words[0]=='A':
            A_list=[float(c) for c in words[1:]]
        elif words[0]=='T':
            CT_list=[float(c) for c in words[1:]]
        elif words[0]=='A_noAwn':
            A_noAwn_list=[float(c) for c in words[1:]]
        elif words[0]=='T_noAwn':
            CT_noAwn_list=[float(c) for c in words[1:]]
        elif words[0]=='spikeLen':
            spikeLen=float(words[1])
        elif words[0]=='spikeletNum':
            spikeletNum=int(words[1])
    
    line=open(fn_weather,'r').readline()
    words=line.strip().split('\t')
    DOY=int(words[1])
    
    fnOut0=op.file_3d[0:-4]+'_ray.txt'
    fnOut0_noAwn=op.file_3d[0:-4]+'_noAwn_ray.txt'
    M_plant_raw_fn=op.file_3d
    M_plant_noAwn_raw_fn=op.file_3d[0:-4]+'_noAwn.txt'
    
    ray_file_list=[]
    ray_file_noAwn_list=[]
    for spikeAng in rotateAng_list:
        M_plant_fn=M_plant_raw_fn[0:-4]+'_'+str(spikeAng)+'.txt'
        M_plant1_fn=M_plant_fn[0:-4]+'_1.txt'
        M_plant2_fn=M_plant_fn[0:-4]+'_2.txt'
        plant3D_1_fn=M_plant_fn[0:-4]+'_1_04-060_045.png'
        M_plant_ray1_fn=M_plant_fn[0:-4]+'_1_ray.txt'
        M_plant_ray2_fn=M_plant_fn[0:-4]+'_2_ray.txt'
        M_plant_ray_fn=M_plant_fn[0:-4]+'_ray.txt'
        
        M_plant_noAwn_fn=M_plant_noAwn_raw_fn[0:-4]+'_'+str(spikeAng)+'.txt'
        M_plant1_noAwn_fn=M_plant_noAwn_fn[0:-4]+'_1.txt'
        M_plant2_noAwn_fn=M_plant_noAwn_fn[0:-4]+'_2.txt'
        plant3D_noAwn_1_fn=M_plant_noAwn_fn[0:-4]+'_1_04-060_045.png'
        M_plant_noAwn_ray1_fn=M_plant_noAwn_fn[0:-4]+'_1_ray.txt'
        M_plant_noAwn_ray2_fn=M_plant_noAwn_fn[0:-4]+'_2_ray.txt'
        M_plant_noAwn_ray_fn=M_plant_noAwn_fn[0:-4]+'_ray.txt'
        
        ray_file_list.append(M_plant_ray_fn)
        ray_file_noAwn_list.append(M_plant_noAwn_ray_fn)
        if not overWrite: 
            print('spikeAng %d Ray tracing (with awn)...'%spikeAng)
            if not os.path.exists(M_plant_ray_fn):
                if not (os.path.exists(M_plant_ray1_fn) and os.path.exists(M_plant_ray2_fn)):
                    if not (os.path.exists(M_plant1_fn) and os.path.exists(M_plant2_fn)):
                        canopyReconstructor(M_plant_raw_fn,img_H,cultivarID,runID,spikeAng,1,fnIn_spikeInfo,overWrite)
                    if not os.path.exists(plant3D_1_fn):
                        plantDrawer(M_plant1_fn,plant3D_1_fn,showOn)
                    PChamberRayTracer(rayDense,M_plant2_fn,M_plant_ray2_fn,M_plant1_fn,M_plant_ray1_fn,fn_weather,DOY)  
                print('Merge ray tracing result of face A and B ...')
                rayFileMerge([M_plant_ray1_fn,M_plant_ray2_fn],M_plant_ray_fn,1)
            
            print('spikeAng %d Ray tracing (without awn)...'%spikeAng)
            if not os.path.exists(M_plant_noAwn_ray_fn):
                if not (os.path.exists(M_plant_noAwn_ray1_fn) and os.path.exists(M_plant_noAwn_ray2_fn)):
                    if not (os.path.exists(M_plant1_noAwn_fn) and os.path.exists(M_plant2_noAwn_fn)):
                        canopyReconstructor(M_plant_noAwn_raw_fn,img_H,cultivarID,runID,spikeAng,0,fnIn_spikeInfo,overWrite)
                    if not os.path.exists(plant3D_noAwn_1_fn):
                        plantDrawer(M_plant1_noAwn_fn,plant3D_noAwn_1_fn,showOn)
                    PChamberRayTracer(rayDense,M_plant2_noAwn_fn,M_plant_noAwn_ray2_fn,M_plant1_noAwn_fn,M_plant_noAwn_ray1_fn,fn_weather,DOY)    
                print('Merge ray tracing result of face A and B ...')
                rayFileMerge([M_plant_noAwn_ray1_fn,M_plant_noAwn_ray2_fn],M_plant_noAwn_ray_fn,1)
        else: # overwrite all possible existing files
            print('spikeAng %d Ray tracing (with awn)...'%spikeAng)
            canopyReconstructor(M_plant_raw_fn,img_H,cultivarID,runID,spikeAng,1,fnIn_spikeInfo,overWrite)
            plantDrawer(M_plant1_fn,plant3D_1_fn,showOn)
            PChamberRayTracer(rayDense,M_plant2_fn,M_plant_ray2_fn,M_plant1_fn,M_plant_ray1_fn,fn_weather,DOY)  
            print('Merge ray tracing result of face A and B ...')
            rayFileMerge([M_plant_ray1_fn,M_plant_ray2_fn],M_plant_ray_fn,1)
            
            print('spikeAng %d Ray tracing (without awn)...'%spikeAng)
            canopyReconstructor(M_plant_noAwn_raw_fn,img_H,cultivarID,runID,spikeAng,0,fnIn_spikeInfo,overWrite)
            plantDrawer(M_plant1_noAwn_fn,plant3D_noAwn_1_fn,showOn)
            PChamberRayTracer(rayDense,M_plant2_noAwn_fn,M_plant_noAwn_ray2_fn,M_plant1_noAwn_fn,M_plant_noAwn_ray1_fn,fn_weather,DOY)    
            print('Merge ray tracing result of face A and B ...')
            rayFileMerge([M_plant_noAwn_ray1_fn,M_plant_noAwn_ray2_fn],M_plant_noAwn_ray_fn,1)
    print('Generating merged ray files: %s and %s'%(fnOut0,fnOut0_noAwn))
    rayData_ave=rayFileMerge(ray_file_list,fnOut0,len(ray_file_list))
    rayData_noAwn_ave=rayFileMerge(ray_file_noAwn_list,fnOut0_noAwn,len(ray_file_noAwn_list))
    rayData_plant_ave=rayData_ave[rayData_ave[:,2]>-6,:] # only keep plant tissues (forbid ground and walls)
    rayData_plant_noAwn_ave=rayData_noAwn_ave[rayData_noAwn_ave[:,2]>-6,:] # only keep plant tissues (forbid ground and walls)

    facet_area1_array=rayData_plant_ave[rayData_plant_ave[:,2]>-5,17] # grain area
    facet_PPFD1_array=rayData_plant_ave[rayData_plant_ave[:,2]>-5,24] # grain PPFD
    facet_area2_array=rayData_plant_ave[rayData_plant_ave[:,2]==-5,17] # awn area
    facet_PPFD2_array=rayData_plant_ave[rayData_plant_ave[:,2]==-5,24] # awn PPFD
    
    facet_area_noAwn_array=rayData_plant_noAwn_ave[:,17] # grain area
    facet_PPFD_noAwn_array=rayData_plant_noAwn_ave[:,24] # grain PPFD
    tot_grain_area=sum(facet_area1_array)
    tot_awn_area=sum(facet_area2_array)
    print('Entire spike surface area: %.2f. De-awned spike surface area: %.2f. Fitting AQ parameters ...'%(tot_grain_area+tot_awn_area,tot_grain_area))
    
    tot_PPF_grain1,tot_PPF_awn1,tot_PPF_grain2,fit_simple,fit_rayTrace=fit_AQ(facet_area1_array,facet_PPFD1_array,facet_area2_array,facet_PPFD2_array,facet_area_noAwn_array,facet_PPFD_noAwn_array,PPFD_list,A_list,CT_list,A_noAwn_list,CT_noAwn_list,Amax_coef,Rd_coef)
    print('Ray tracing fitting result: ')
    print('    theta,Amax,phi_CO2,Rd,theta_uniq,Amax_uniq,phi_CO2_uniq,fit_error,R2_fit: ',fit_rayTrace[0:-1])
    
    fnOut=open(fnOut,'w')
    fnOut.write('Description\ttheta\tAmax_gross25\tAmax_net25\tphi_CO2\tRd25\ttheta_awn\tAmax_gross25_awn\tAmax_net25_awn\tphi_CO2_awn\tRd25_awn\ttheta_uniq\tAmax_uniq\tphi_CO2_uniq\tRd_uniq\ttheta_awn_uniq\tAmax_awn_uniq\tphi_CO2_awn_uniq\tRd_awn_uniq\tfit_error\tR2_fit\n')
    content=[fit_rayTrace[0]]+[fit_rayTrace[1]+fit_rayTrace[3]]+fit_rayTrace[1:4]+[fit_rayTrace[4]]+[fit_rayTrace[5]+fit_rayTrace[7]]+fit_rayTrace[5:-1]
    fnOut.write('rayTracing_fit\t'+'\t'.join([str(c) for c in content])+'\n')
    fnOut.write('Surface_area(cm2)\t'+str(round(tot_grain_area+tot_awn_area,2))+'\t'+str(round(tot_grain_area,2))+'\n')
    content=[tot_PPF_grain1[0]/(tot_grain_area+tot_awn_area)*10,tot_PPF_awn1[0]/tot_awn_area*10,tot_PPF_grain2[0]/tot_grain_area*10]
    fnOut.write('PPFD_ave(umol/m2/s)\t'+'\t'.join([str(round(c,0)) for c in content])+'\n')
    fnOut.write('PPFD_in(umol/m2/s)\t'+'\t'.join([str(c) for c in PPFD_list+PPFD_list])+'\n')
    fnOut.write('A_true(nmol/s)\t'+'\t'.join([str(c) for c in A_list+A_noAwn_list])+'\n')
    fnOut.write('tot_PPF_grain(nmol/s)\t'+'\t'.join([str(c) for c in np.hstack((tot_PPF_grain1,tot_PPF_grain2))])+'\n')
    fnOut.write('tot_PPF_awn(nmol/s)\t'+'\t'.join([str(c) for c in tot_PPF_awn1])+'\n')
    fnOut.write('rayTracing_fitA\t'+'\t'.join([str(round(c,2)) for c in fit_rayTrace[-1]])+'\n')
    
    fnOut.close()
    eclipse=time.time()-start
    print('spikeAQ fitting done! Time used: %.1f s.'%eclipse)