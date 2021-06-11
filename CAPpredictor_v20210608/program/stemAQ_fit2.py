#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
stemAQ_fit2.py

1. construct 3D wheat stem
2. perform ray tracing simulation
3. fit AQ curve
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

# interpolant = PchipInterpolator(x, y)
# xnew = np.linspace(0, 10, num=41, endpoint=True)
# interpolant(xnew)

__author__="CHANG Tiangen"
__date__="20201215"
__copyright__="Copyright 2020 Aspar CHANG"
__license__="SIPPE"
__version__="0.1"
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
    description_string='The program is going to fit the wheat AQ curve.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--fnIn', dest='file_in', help='Input stem basic info file')
    parser.add_argument('-iw', '--fnIn_w', dest='file_in_w', help='Input weather data file')
    parser.add_argument('-o1', '--fnOut1', dest='file_out1', help='Output structure file')
    parser.add_argument('-o2', '--fnOut2', dest='file_out2', help='Output AQ parameter file')
    parser.add_argument('-rID', '--runID', dest='runID', type=int, help='Remark of current run ID, used as file prefix and random seed.')
    parser.add_argument('-log', '--printInfo2logFile', dest='log', type=int, default=0, help='Log printed info to log file (1) or not (0).')
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
    #print(popt_new,popt_old)
    theta,Amax25,phi_CO2,Rd25=popt_old
    theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq=uniq_old
    if fit_error_new<fit_error-epsilon:
        print('fit_error_new: %.2f'%fit_error_new)
        popt_new=popt
        fit_error_out=fit_error_new
        itera_test_out=1
        theta_uniq='YES'
        Amax25_uniq='YES'
        phi_CO2_uniq='YES'
        Rd25_uniq='YES'
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
    uniq_new=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq]
    return itera_test_out,fit_error_out,popt_new,uniq_new

def fit_AQ(facet_area_array,facet_PPFD_array,PAR_ave,A_ave,CT_ave,Amax_coef,Rd_coef):
    x=np.array(PAR_ave)
    y=np.array(A_ave)
    z=np.array(CT_ave)
    tot_PPF=np.dot(facet_PPFD_array,facet_area_array)/10 # nmol s-1
    x_tot=np.array(PAR_ave)/2000*tot_PPF
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
    
    while itera_test < 50:
        print('            Simple AQ fitting: iteration %d, fit_error %.1f'%(itera_test,fit_error))
        theta_ini=random.uniform(0.01,0.999)
        Amax25_ini=random.randint(1,100)*1.01
        phi_CO2_ini=random.uniform(0,0.1)
        Rd25_ini=random.randint(1,100)*1.01
        paras0=np.array([theta_ini,Amax25_ini,phi_CO2_ini,Rd25_ini]) # theta,Amax,alpha,Rd
        popt= fmin_slsqp(fit_residuals, paras0, bounds=[(0.01,0.999),(0,101),(0,0.1),(0,101)], args=(x_tot,y,z,Amax_coef,Rd_coef),iprint = 0)
        popt=popt[0:]
        if math.isnan(popt[0]):
            continue
        if theta+Amax25+phi_CO2==0:
            theta,Amax25,phi_CO2,Rd25=popt
            fit_error=fit_residuals(popt[0:],x_tot,y,z,Amax_coef,Rd_coef)
        else:
            fit_error_new=fit_residuals(popt[0:],x_tot,y,z,Amax_coef,Rd_coef)
            popt_old=[theta,Amax25,phi_CO2,Rd25]
            uniq_old=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq]
            itera_test,fit_error,popt_new,uniq_new = uniqueFitJudger(itera_test,fit_error_new,fit_error,popt,popt_old,uniq_old)
            theta,Amax25,phi_CO2,Rd25=popt_new
            theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq=uniq_new
        itera_test=itera_test+1
    Rd_array=Rd25*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Amax_array=Amax25*(Amax_coef[0]*z**2 + Amax_coef[1]*z + Amax_coef[2]) # Rd-T empirical relation
    Ag_max=Amax_array+Rd_array
    A_fit_array=(phi_CO2*x_tot+Ag_max-np.sqrt((phi_CO2*x_tot+Ag_max)**2-4*theta*phi_CO2*x_tot*Ag_max))/(2*theta)-Rd_array
    
    R2_fit=1-sum((A_fit_array-y)**2)/sum((y-np.mean(y))**2)
    result_1=[theta,Amax25,phi_CO2,Rd25,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,fit_error,R2_fit,A_fit_array]
    
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
    Amax25_ub=20
    Rd25_ub=10
    while itera_test < 50:
        print('            Ray tracing AQ fitting: iteration %d, fit_error %.1f'%(itera_test,fit_error))
        theta_ini=random.uniform(0.01,0.999)
        Amax25_ini=random.uniform(0,Amax25_ub)
        phi_CO2_ini=random.uniform(0,0.1)
        Rd25_ini=random.uniform(0,Rd25_ub)
#         theta_ini=0.839100165
#         Amax_ini=4.093976743
#         phi_CO2_ini=0.079715078
#         Rd25_ini=2.420171864
        paras0=np.array([theta_ini,Amax25_ini,phi_CO2_ini,Rd25_ini]) # theta,Amax,alpha,Rd
        popt= fmin_slsqp(fit_residuals2, paras0, bounds=[(0.01,0.999),(0,Amax25_ub),(0,0.1),(0,Rd25_ub)], args=(facet_area_array,facet_PPFD_array,x,y,z,Amax_coef,Rd_coef),iprint = 0)
        popt=popt[0:]
        if math.isnan(popt[0]):
            continue
        if theta+Amax25+phi_CO2==0:
            theta,Amax25,phi_CO2,Rd25=popt
            fit_error=fit_residuals2(popt[0:],facet_area_array,facet_PPFD_array,x,y,z,Amax_coef,Rd_coef)
        else:
            fit_error_new=fit_residuals2(popt[0:],facet_area_array,facet_PPFD_array,x,y,z,Amax_coef,Rd_coef)
            popt_old=[theta,Amax25,phi_CO2,Rd25]
            uniq_old=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq]
            itera_test,fit_error,popt_new,uniq_new = uniqueFitJudger(itera_test,fit_error_new,fit_error,popt,popt_old,uniq_old)
            theta,Amax25,phi_CO2,Rd25=popt_new
            theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq=uniq_new
        itera_test=itera_test+1
    AQ_paras=[theta,Amax25,phi_CO2,Rd25]
    A_fit_array2=Aspike_fit(facet_area_array,facet_PPFD_array,AQ_paras,x,z,Amax_coef,Rd_coef)
    R2_fit=1-sum((A_fit_array2-y)**2)/sum((y-np.mean(y))**2)
    result_2=[theta,Amax25,phi_CO2,Rd25,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,fit_error,R2_fit,A_fit_array2]

    return x_tot,result_1,result_2

def Aspike_fit(facet_area_array,facet_PPFD_array,AQ_paras,PAR_ave,CT_ave,Amax_coef,Rd_coef):
    theta,Amax25,phi_CO2,Rd25=AQ_paras
    y_fit=[]
    for PAR_i in range(len(PAR_ave)):
        x=facet_PPFD_array*PAR_ave[PAR_i]/2000
        Rd=Rd25*(Rd_coef[0]*CT_ave[PAR_i]**2 + Rd_coef[1]*CT_ave[PAR_i] + Rd_coef[2]) # Rd-T empirical relation
        Amax=Amax25*(Amax_coef[0]*CT_ave[PAR_i]**2 + Amax_coef[1]*CT_ave[PAR_i] + Amax_coef[2]) # Amax-T empirical relation
        Ag_max=Amax+Rd
        A_fit_facet_array=(phi_CO2*x+Ag_max-np.sqrt((phi_CO2*x+Ag_max)**2-4*theta*phi_CO2*x*Ag_max))/(2*theta)-Rd
        y_fit.append(float(np.dot(facet_area_array,A_fit_facet_array)/10)) # nmol ear-1 s-1
    return np.array(y_fit)

def fit_residuals(paras,x_tot,y,z,Amax_coef,Rd_coef):#y=A_entire, x=PPFD, z=CT, y2=A_noAwn, z2=CT2
    theta,Amax25,phi_CO2,Rd25=paras
    Rd_facet_array=Rd25*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Amax_facet_array=Amax25*(Amax_coef[0]*z**2 + Amax_coef[1]*z + Amax_coef[2]) # A-T empirical relation
    Ag_facet_max=Amax_facet_array+Rd_facet_array
    
    sqrt1=(phi_CO2*x_tot+Ag_facet_max)**2-4*theta*phi_CO2*x_tot*Ag_facet_max
    y_fit=(phi_CO2*x_tot+Ag_facet_max-np.sqrt(sqrt1))/(2*theta)-Rd_facet_array
    fit_err=sum((y_fit-y)**2)
    #print('Simple fit: ',y_fit,y_entire)
    return fit_err
    
def fit_residuals2(paras,facet_area_array,facet_PPFD_array,PAR_ave,A_ave,CT_ave,Amax_coef,Rd_coef):
    theta,Amax25,phi_CO2,Rd25=paras
    Rd_facet_array=Rd25*(Rd_coef[0]*CT_ave**2 + Rd_coef[1]*CT_ave + Rd_coef[2]) # Rd-T empirical relation    
    Amax_facet_array=Amax25*(Amax_coef[0]*CT_ave**2 + Amax_coef[1]*CT_ave + Amax_coef[2]) # A-T empirical relation
    y_fit=[]
    for PAR_i in range(len(PAR_ave)):
        x_facet=facet_PPFD_array*(PAR_ave[PAR_i]/2000)
        Rd=Rd_facet_array[PAR_i]
        Amax=Amax_facet_array[PAR_i]
        Ag_max=Amax+Rd        
        sqrt1=(phi_CO2*x_facet+Ag_max)**2-(4*theta*phi_CO2*Ag_max)*x_facet
        A_fit_array=(phi_CO2*x_facet+Ag_max-np.sqrt(sqrt1))/(2*theta)-Rd
        y_fit.append(float(np.dot(facet_area_array,A_fit_array)/10)) # nmol ear-1 s-1
    #print('rayTrace fit: ',y_fit,A_entire)
    fit_err=sum((np.array(y_fit)-A_ave)**2)
    return fit_err

def stemReconstructor(M_plant_raw_fn,runID,stemLen,stemDiameter,stemN,overWrite):
    # Canopy structure reconstructing
    print('    Stem structure reconstructing ...\n')
    command='python .\pchamberStemConstructor.py -o %s -rID %d -ow %d -log 1 -sl %.2f -sd %.3f -n %.2f'%(M_plant_raw_fn,runID,overWrite,stemLen,stemDiameter,stemN)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command) 
    
def plantDrawer(M_plant_fn,plant3D_fn,showOn):
    print('    Stem structure plotting ...\n')
    command='python .\plantDraw3D.py -i %s -o %s -ccs 4 -cql 0 -cqu 1 -cm plasma -alpha 1 -va1 60 -show %d'%(M_plant_fn,plant3D_fn,showOn)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command) 

def PChamberRayTracer(rayDense,M_plant2_fn,M_plant_ray2_fn,M_plant1_fn,M_plant_ray1_fn,fn_weather,DOY):
    print('P-Chamber stem ray tracing. Face A ...\n')
    command='.\\fastTracerVS2019_0.exe -L 31 -d %d -D -0.01 29.01 -2.01 2.01 0 5.02 -S 12 -n %f -W 6 2 7 -m %s -o %s -C %s -s 1'%(DOY,rayDense,M_plant1_fn,M_plant_ray1_fn,fn_weather)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)   
    print('P-Chamber stem ray tracing. Face B ...\n')
    command='.\\fastTracerVS2019_0.exe -L 31 -d %d -D -0.01 29.01 -2.01 2.01 0 5.02 -S 12 -n %f -W 6 2 7 -m %s -o %s -C %s -s 1'%(DOY,rayDense,M_plant2_fn,M_plant_ray2_fn,fn_weather)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)   

def file_rename(fnOut_str,runID):
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
    return fnOut

if __name__=="__main__":
    start=time.time()
    op=createHelp()
    log2file=op.log
    if log2file:
        sys.stdout = Logger('./log_stemAQ_fit')
    width_format=94
    welcome_info='Welcome to use stem AQ fitter.'
    time_info=time.asctime( time.localtime(time.time()) )
    print(welcome_info.center(width_format, '-'))
    print(time_info.center(width_format, ' '))
    
    runID=op.runID
    random.seed(runID)
    fn_weather=op.file_in_w
    fnIn_spikeInfo=op.file_in
    fnOut_str1=op.file_out1
    fnOut_str2=op.file_out2
    overWrite=op.ow
    rayDense=op.RD
    showOn=op.showOn
    
    fnOut1=file_rename(fnOut_str1,runID)
    fnOut=file_rename(fnOut_str2,runID)
    
    op.file_out1=fnOut1
    op.file_out2=fnOut
    data = json.dumps(vars(op), indent=4,ensure_ascii=False, sort_keys=False,separators=(',', ':'))
    print('Running parameters:')
    print(data)
    
    Rd_coef=[0,0,1]
    Amax_coef=[0,0,1]
    CT_list=[25]
    for line in open(fnIn_spikeInfo,'r').readlines():
        words=line.strip().split('\t')
        if words[0]=='PPFD':
            PPFD_list=[int(c) for c in words[1:]]
        elif words[0]=='A':
            A_list=[float(c) for c in words[1:]] # morning
        elif words[0]=='CT':
            CT_list=[float(c) for c in words[1:]] # morning
        elif words[0]=='stemLen':
            stemLen=float(words[1])
        elif words[0]=='stemDiameter':
            stemDiameter=float(words[1])*1.047 # *1.047 to ensure Perimeter of hexagon equals to circle (area ratio 91%)
        elif words[0]=='stemN':
            stemN=float(words[1])
        elif words[0]=='coef_stem_Asat_airT':
            Amax_coef=[float(c) for c in words[1:]]
        elif words[0]=='coef_stem_Rd_airT':
            Rd_coef=[float(c) for c in words[1:]]
    if len(CT_list)==1:
        CT_list=CT_list*len(PPFD_list)
    if (len(PPFD_list)!=len(A_list)) or (len(PPFD_list)!=len(CT_list)):
        raise Exception('ERROR: length of PPFD_list, A_list and CT_list not match in file: %s'%fnIn_spikeInfo)
    
    line=open(fn_weather,'r').readline()
    words=line.strip().split('\t')
    DOY=int(words[1])
    
    fnOut0=fnOut1[0:-4]+'_ray.txt'
    M_plant_raw_fn=fnOut1
    
    ray_file_list=[]
    M_plant_fn=M_plant_raw_fn[0:-4]+'.txt'
    M_plant1_fn=M_plant_fn[0:-4]+'_1.txt'
    M_plant2_fn=M_plant_fn[0:-4]+'_2.txt'
    plant3D_1_fn=M_plant_fn[0:-4]+'_1_04-060_045.png'
    M_plant_ray1_fn=M_plant_fn[0:-4]+'_1_ray.txt'
    M_plant_ray2_fn=M_plant_fn[0:-4]+'_2_ray.txt'
    M_plant_ray_fn=M_plant_fn[0:-4]+'_ray.txt'
    
    ray_file_list.append(M_plant_ray_fn)
    if not overWrite: 
        print('Ray tracing ...')
        if not os.path.exists(M_plant_ray_fn):
            if not (os.path.exists(M_plant_ray1_fn) and os.path.exists(M_plant_ray2_fn)):
                if not (os.path.exists(M_plant1_fn) and os.path.exists(M_plant2_fn)):
                    stemReconstructor(M_plant_raw_fn,runID,stemLen,stemDiameter,stemN,overWrite)
                if not os.path.exists(plant3D_1_fn):
                    plantDrawer(M_plant1_fn,plant3D_1_fn,showOn)
                PChamberRayTracer(rayDense,M_plant2_fn,M_plant_ray2_fn,M_plant1_fn,M_plant_ray1_fn,fn_weather,DOY)  
            print('Merge ray tracing result of face A and B ...')
            rayFileMerge([M_plant_ray1_fn,M_plant_ray2_fn],M_plant_ray_fn,1)
    else: # overwrite all possible existing files
        print('Ray tracing ...')
        stemReconstructor(M_plant_raw_fn,runID,stemLen,stemDiameter,stemN,overWrite)
        plantDrawer(M_plant1_fn,plant3D_1_fn,showOn)
        PChamberRayTracer(rayDense,M_plant2_fn,M_plant_ray2_fn,M_plant1_fn,M_plant_ray1_fn,fn_weather,DOY)  
        rayFileMerge([M_plant_ray1_fn,M_plant_ray2_fn],M_plant_ray_fn,1)
    print('Generating merged ray files: %s'%(fnOut0))
    rayData_ave=rayFileMerge(ray_file_list,fnOut0,len(ray_file_list))
    rayData_plant_ave=rayData_ave[rayData_ave[:,2]>-6,:] # only keep plant tissues (forbid ground and walls)
    facet_area1_array=rayData_plant_ave[:,17] # grain area
    facet_PPFD1_array=rayData_plant_ave[:,24] # grain PPFD
    print('Fitting AQ parameters ...')
    #raise Exception(A_list,CT_list,Amax_coef,Rd_coef)
    tot_PPF_stem,fit_simple,fit_rayTrace=fit_AQ(facet_area1_array,facet_PPFD1_array,PPFD_list,A_list,CT_list,Amax_coef,Rd_coef)
    tot_area=sum(facet_area1_array)
    mean_PPFD_max=tot_PPF_stem[0]/tot_area*10
    print('Stem surface area (cm2): %.2f. '%tot_area)
    print('Stem average PPFD (umol/m2/s): %d. '%mean_PPFD_max)
    print('Ray tracing fit result: ')
    print('    theta,Amax,phi_CO2,Rd,theta_uniq,Amax_uniq,phi_CO2_uniq,fit_error,R2_fit: ',fit_rayTrace[0:-1])
    
    fnOut=open(fnOut,'w')
    fnOut.write('Description\ttheta\tAmax25\tAmax25_net\tphi_CO2\tRd25\ttheta_uniq\tAmax25_uniq\tphi_CO2_uniq\tRd_uniq\tfit_error\tR2_fit\n')
    content=[fit_rayTrace[0]]+[fit_rayTrace[1]+fit_rayTrace[3]]+fit_rayTrace[1:-1]
    fnOut.write('rayTracing_fit\t'+'\t'.join([str(c) for c in content])+'\n')
    fnOut.write('Surface_area(cm2)\t'+str(round(tot_area,2))+'\n')
    fnOut.write('PPFD_ave(umol/m2/s)\t'+str(round(mean_PPFD_max,0))+'\n')
    fnOut.write('PPFD_in(umol/m2/s)\t'+'\t'.join([str(c) for c in PPFD_list])+'\n')
    fnOut.write('A_true(nmol/s)\t'+'\t'.join([str(c) for c in A_list])+'\n')
    fnOut.write('tot_PPF_stem(nmol/s)\t'+'\t'.join([str(c) for c in tot_PPF_stem])+'\n')
    fnOut.write('rayTracing_fitA\t'+'\t'.join([str(round(c,2)) for c in fit_rayTrace[-1]])+'\n')
    
    fnOut.close()
    eclipse=time.time()-start
    print('stemAQ fitting done! Time used: %.1f s.'%eclipse)