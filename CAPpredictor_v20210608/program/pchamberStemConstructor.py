#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
pchamberStemConstructor.py

construct 3D stem in p-chamber for ray tracing algorithm
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

import matplotlib.pyplot as plt
#from scipy.interpolate import pchip
from scipy.interpolate import PchipInterpolator

# interpolant = PchipInterpolator(x, y)
# xnew = np.linspace(0, 10, num=41, endpoint=True)
# interpolant(xnew)

__author__="CHANG Tiangen"
__date__="20201216"
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
    description_string='The program is going to construct a 3D stem in P-chamber.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-o', '--output file', dest='file_out', help='Output file')
    parser.add_argument('-rID', '--run ID', dest='runID', type=int, help='Remark of current run ID, used as file prefix and random seed.')
    parser.add_argument('-log', '--printInfo2logFile', dest='log', type=int, default=1, help='Log printed info to log file (1) or not (0).')
    parser.add_argument('-ow', '--overWrite', dest='ow', type=int, default=0, help='Overwrite the existing files (1) or not (0).')
    parser.add_argument('-sl', '--stemLen', dest='stemLen', type=float, help='Stem length (cm).')
    parser.add_argument('-sd', '--stemDiameter', dest='stemDiameter', type=float, help='Stem Diameter (cm).')
    parser.add_argument('-n', '--N', dest='stemN', type=float, help='Stem nitrogen content (g/g).')

    op=parser.parse_args()
    return op

def buildStem(barcode,LNC_Kt_Kr,organ_base,organ_tip,organ_diameter):
    x_in_order=[1,1/2,-1/2,-1,-1/2,1/2,1]
    y_in_order=[0,math.sqrt(3)*1/2,math.sqrt(3)*1/2,0,-math.sqrt(3)*1/2,-math.sqrt(3)*1/2,0]
    z_vec=np.array([0,0,1])
    M_stem_all=np.zeros([0,17])
    organ_vec=np.array(organ_tip)-organ_base
    organ_len=np.linalg.norm(organ_vec)
    seg_number=int(round(organ_len))+1 # cut stem into pieces with len < 1cm
    seg_len=organ_len/seg_number    
    organ_tip_temp=[0,0,0]
    for pt_i in range(seg_number):
        organ_base_temp=organ_tip_temp
        organ_tip_temp=organ_base_temp+z_vec*seg_len
        x0 = organ_base_temp[0]
        y0 = organ_base_temp[1]
        z0 = organ_base_temp[2]
        dem0 = organ_diameter[0]*(organ_len-seg_len*pt_i)/organ_len+organ_diameter[1]*seg_len*pt_i/organ_len
        radi0 = dem0/2
        x1 = organ_tip_temp[0]
        y1 = organ_tip_temp[1]
        z1 = organ_tip_temp[2]
        dem1 = organ_diameter[0]*(organ_len-seg_len*(pt_i+1))/organ_len+organ_diameter[1]*seg_len*(pt_i+1)/organ_len
        radi1 = dem1/2
        for i in range(6):
            p_s_01 = [radi0*x_in_order[i] +   x0, radi0*y_in_order[i] +   y0, z0]
            p_s_02 = [radi0*x_in_order[i+1] + x0, radi0*y_in_order[i+1] + y0, z0]
            p_s_11 = [radi1*x_in_order[i] +   x1, radi1*y_in_order[i] +   y1, z1]
            p_s_12 = [radi1*x_in_order[i+1] + x1, radi1*y_in_order[i+1] + y1, z1]
            M_stem = [barcode + p_s_01 + p_s_02 + p_s_11 + LNC_Kt_Kr] + [barcode + p_s_02 + p_s_11 + p_s_12 + LNC_Kt_Kr]
            #print(M_stem)
            M_stem_all = np.vstack((M_stem_all,M_stem))
    organ_direction=(organ_vec)/organ_len
    rotation_vec=np.cross(z_vec,organ_direction)
    sin_a=np.linalg.norm(rotation_vec)
    triangle_num=np.size(M_stem_all,0)
    off_set=np.tile(organ_base,(triangle_num,1))
    if abs(sin_a)>0.001:
        cos_a=np.dot(z_vec,organ_direction)
        rotation_vec=rotation_vec/sin_a
        a=rotation_vec[0]
        b=rotation_vec[1]
        c=rotation_vec[2]
        row_1=[a**2+(1-a**2)*cos_a,a*b*(1-cos_a)+c*sin_a,a*c*(1-cos_a)-b*sin_a]
        row_2=[a*b*(1-cos_a)-c*sin_a,b**2+(1-b**2)*cos_a,b*c*(1-cos_a)+a*sin_a]
        row_3=[a*c*(1-cos_a)+b*sin_a,c*b*(1-cos_a)-a*sin_a,c**2+(1-c**2)*cos_a]
        r_mat=np.array([row_1,row_2,row_3])
    else:
        r_mat=np.eye(3)
    M_stem_all[:,5:8]=off_set + np.dot(M_stem_all[:,5:8],r_mat)
    M_stem_all[:,8:11]=off_set + np.dot(M_stem_all[:,8:11],r_mat)
    M_stem_all[:,11:14]=off_set + np.dot(M_stem_all[:,11:14],r_mat)
    return M_stem_all

def buildChamber(M_plant,Kr_wall,chamber_H,chamber_W,chamber_L):
    cos_xRotAng=-1
    sin_xRotAng=0
    x_rot_mat=np.array([[1,0,0],[0,cos_xRotAng,-sin_xRotAng],[0,sin_xRotAng,cos_xRotAng]])
    r_mat2=x_rot_mat
    
    z_move=[0,0,chamber_H/2]
    triangle_num=np.size(M_plant,0)
    off_set=np.tile(z_move,(triangle_num,1))
    M_tiller=np.zeros([triangle_num,17])
    M_tiller[:,0:5]=M_plant[:,0:5]
    M_tiller[:,14:]=M_plant[:,14:]
    M_tiller[:,5:8]=off_set + M_plant[:,5:8]
    M_tiller[:,8:11]=off_set + M_plant[:,8:11]
    M_tiller[:,11:14]=off_set + M_plant[:,11:14]
    
    M_chamber1=np.zeros([triangle_num+1000,17])
    M_chamber2=np.zeros([triangle_num+1000,17])
    current_pos=triangle_num
    M_chamber1[0:current_pos,:]=M_tiller
    
    M_tiller=np.zeros([triangle_num,17])
    M_tiller[:,0:5]=M_plant[:,0:5]
    M_tiller[:,14:]=M_plant[:,14:]
    M_tiller[:,5:8]=off_set + np.dot(M_plant[:,5:8],r_mat2)
    M_tiller[:,8:11]=off_set + np.dot(M_plant[:,8:11],r_mat2)
    M_tiller[:,11:14]=off_set + np.dot(M_plant[:,11:14],r_mat2)
    M_chamber2[0:current_pos,:]=M_tiller
    # z move chamber_H/2
    xmin=0
    ymin=-chamber_W/2
    zmin=0.01
    xmax=chamber_L
    ymax=chamber_W/2
    zmax=chamber_H+zmin
    
    x_steps=np.linspace(xmin,xmax,round(xmax-xmin)+1)[0:-1]
    y_steps=np.linspace(ymin,ymax,round(ymax-ymin)+1)[0:-1]
    z_steps=np.linspace(zmin,zmax,round(zmax-zmin)+1)[0:-1]
    # construct the P-chamber wall
    for x_i in x_steps:
        for y_i in y_steps:
            p1=[x_i,y_i,zmin]
            p2=[x_i,y_i+1,zmin]
            p3=[x_i+1,y_i,zmin]
            p4=[x_i+1,y_i+1,zmin]
            M_wall = np.array([[0,0,-8,1,0] + p1 + p2 + p4 + [0,0,Kr_wall]] + [[0,0,-8,1,0] + p1 + p4 + p3 + [0,0,Kr_wall]])
            M_chamber1[current_pos:current_pos+2,:]=M_wall
            current_pos+=2
    for z_i in z_steps:
        for y_i in y_steps:
            p1=[xmin,y_i,z_i]
            p2=[xmin,y_i,z_i+1]
            p3=[xmin,y_i+1,z_i]
            p4=[xmin,y_i+1,z_i+1]
            M_wall = np.array([[0,0,-8,1,0] + p1 + p2 + p4 + [0,0,Kr_wall]] + [[0,0,-8,1,0] + p1 + p4 + p3 + [0,0,Kr_wall]])
            M_chamber1[current_pos:current_pos+2,:]=M_wall
            current_pos+=2
    for z_i in z_steps:
        for y_i in y_steps:
            p1=[xmax,y_i,z_i]
            p2=[xmax,y_i,z_i+1]
            p3=[xmax,y_i+1,z_i]
            p4=[xmax,y_i+1,z_i+1]
            M_wall = np.array([[0,0,-8,1,0] + p1 + p4 + p2 + [0,0,Kr_wall]] + [[0,0,-8,1,0] + p1 + p3 + p4 + [0,0,Kr_wall]])
            M_chamber1[current_pos:current_pos+2,:]=M_wall
            current_pos+=2
    for z_i in z_steps:
        for x_i in x_steps:
            p1=[x_i,ymin,z_i]
            p2=[x_i,ymin,z_i+1]
            p3=[x_i+1,ymin,z_i]
            p4=[x_i+1,ymin,z_i+1]
            M_wall = np.array([[0,0,-8,1,0] + p1 + p4 + p2 + [0,0,Kr_wall]] + [[0,0,-8,1,0] + p1 + p3 + p4 + [0,0,Kr_wall]])
            M_chamber1[current_pos:current_pos+2,:]=M_wall
            current_pos+=2
    for z_i in z_steps:
        for x_i in x_steps:
            p1=[x_i,ymax,z_i]
            p2=[x_i,ymax,z_i+1]
            p3=[x_i+1,ymax,z_i]
            p4=[x_i+1,ymax,z_i+1]
            M_wall = np.array([[0,0,-8,1,0] + p1 + p2 + p4 + [0,0,Kr_wall]] + [[0,0,-8,1,0] + p1 + p4 + p3 + [0,0,Kr_wall]])
            M_chamber1[current_pos:current_pos+2,:]=M_wall
            current_pos+=2
    M_chamber2[triangle_num:current_pos,:]=M_chamber1[triangle_num:current_pos,:]
    return M_chamber1[0:current_pos,:],M_chamber2[0:current_pos,:]

def SLN2Chl(SLN):
    # gN/m2. -> ug/cm2. this equation combined two high-cited papers, the logic is SPAD -> LN (paper1); SPAD -> Chl (paper2); then LN -> Chl (solve eqn.)
    try:
        Chl=7.252*SLN**2+18.191*SLN+1.5012
    except:
        Chl=list(map(lambda x:7.252*x**2+18.191*x+1.5012,SLN))
    return Chl

def Chl2Kr(Chl):
    # fitting from PROSPECT5 model. http://opticleaf.ipgp.fr/index.php?page=prospect. 
    try:
        Kr=0.3605*Chl**(-0.502)
    except:
        Kr=list(map(lambda x:0.3605*x**(-0.502),Chl))
    return Kr

def Chl2Kt(Chl):
    # fitting from PROSPECT5 model. http://opticleaf.ipgp.fr/index.php?page=prospect. 
    try:
        Kt=max(0,-0.082*math.log(Chl) + 0.3761)*1
    except:
        Kt=list(map(lambda x:max(0,-0.082*math.log(x) + 0.3761),Chl))
    return Kt

    
if __name__=="__main__":
    start=time.time()
    
    op=createHelp()
    log2file=op.log
    if log2file:
        sys.stdout = Logger('./log_pchamberStemConstructor')
    width_format=94
    welcome_info='Welcome to use wheat canopy constructor'
    time_info=time.asctime( time.localtime(time.time()) )
    print(welcome_info.center(width_format, '-'))
    print(time_info.center(width_format, ' '))
    
    runID=op.runID
    random.seed(runID)
    fnOut_str=op.file_out
    overWrite=op.ow
    stemLen=op.stemLen
    stemDiameter=op.stemDiameter
    stemN=op.stemN
    fnOut_list=fnOut_str.split('\\')
    resultDir='\\'.join(fnOut_list[0:-1])
    if not resultDir:
        resultDir='.\\'
    elif not os.path.exists(resultDir):
        print('  Result directory not exists, creating it ...')
        try:
            os.mkdir(resultDir)
        except:
            raise Exception('Cannot create result directory: %s'%resultDir)
    if 'Run' not in fnOut_list[-1]:
        fnOut_list[-1]='Run'+str(runID)+'_'+fnOut_list[-1]
    fnOut='\\'.join(fnOut_list)
    op.file_out=fnOut
    
    data = json.dumps(vars(op), indent=4,ensure_ascii=False, sort_keys=False,separators=(',', ':'))
    print('Running parameters:')
    print(data)
    
    fnOut_pchamber1=fnOut[0:-4]+'_1.txt'
    fnOut_pchamber2=fnOut[0:-4]+'_2.txt'
    
    if (not overWrite) and os.path.exists(fnOut_pchamber1)  and os.path.exists(fnOut_pchamber2):
        print('  Existing 3D canopy structure file "%s" found! Nothing more needs to do, exit directly.'%fnOut)
        sys.exit()
    
    stem_LNC=stemN
    stem_SLN=stem_LNC/100*43 # gN/g * g/m2
    stem_Chl=SLN2Chl(stem_SLN)
    Kr_stem=Chl2Kr(stem_Chl)/2 # if not, when upper stem =1.05 g/g N -> Kr=0.107, too high
    Kt_stem=0
    
    barcode=[1,1,0,3,0]
    LNC_Kt_Kr=[stem_LNC,Kt_stem,Kr_stem]
    organ_base=[0,0,0]
    organ_tip=[stemLen,0,0]
    M_stem = buildStem(barcode,LNC_Kt_Kr,organ_base,organ_tip,[stemDiameter,stemDiameter]) # plant, tiller, Position, organCode
    print('********Stem length, stem diameter: %.2f cm, %.3f cm. ********'%(stemLen,stemDiameter))
    # Spike in P Chamber   
    chamber_H=5
    chamber_W=4
    chamber_L=29
    Kr_wall=0.75
    [M_canopy1,M_canopy2] = buildChamber(M_stem,Kr_wall,chamber_H,chamber_W,chamber_L)
    print('  Writing Triangle-Patch into file...')
    format_str='%d\t'*5+'%.3f\t'*12
    format_str=format_str[0:-1]
    np.savetxt(fnOut_pchamber1,M_canopy1,fmt=format_str) # , delimiter='\t'
    np.savetxt(fnOut_pchamber2,M_canopy2,fmt=format_str) # , delimiter='\t'
    eclipse=time.time()-start
    print('  Stem-in-chamber construction done! Time used: ', eclipse)
    