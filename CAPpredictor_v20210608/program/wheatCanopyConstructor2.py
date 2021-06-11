#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
wheatCanopyConstructor2.py

construct 3D wheat canopy with 2D photos for ray tracing algorithm
Update to 2: Kr and Kt are determined by: LNC -> SPAD -> chl
Remark for this script:
   spike internode length = 0
   leaves are straightened
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
# interpolant = PchipInterpolator(x, y)
# xnew = np.linspace(0, 10, num=41, endpoint=True)
# interpolant(xnew)

__author__="CHANG Tiangen"
__date__="20201230"
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
    description_string='The program is going to construct a 3D wheat canopy.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-di', '--input directory', dest='dir_input', help='Directory for input files')
    parser.add_argument('-i', '--input file', dest='file_in', help='Input raw data file of manually-extracted tiller structure')
    parser.add_argument('-i2', '--input file2', dest='file_in2', nargs='+',help='Input configure files (separated by space)')
    parser.add_argument('-o', '--output file', dest='file_out', help='Output file')
    parser.add_argument('-ih', '--image height', dest='img_H', type=int, help='Image height in pixels')
    parser.add_argument('-lineNA', '--cultivar line name', dest='lineNA', type=str, help='Cultivar specific key word in image file name (N or Y)')
    parser.add_argument('-sp', '--spike', dest='withSpike', type=int,default=1, help='Canopy with spike (1) or de-spiked (0)')
    parser.add_argument('-rID', '--run ID', dest='runID', type=int, help='Remark of current run ID, used as file prefix and random seed.')
    parser.add_argument('-log', '--printInfo2logFile', dest='log', type=int, default=0, help='Log printed info to log file (1) or not (0).')
    parser.add_argument('-ow', '--overWrite', dest='ow', type=int, default=0, help='Overwrite the existing files (1) or not (0).')

    op=parser.parse_args()
    return op

def inputFileProcessing(dataIn,img_H,keyWord,spikeInternodeLenZero,leafStraighten):
    print('  Input key point marking file processing...')
    dataIn_len=len(dataIn)
    scale=1
    i=0
    plant_info_str=''
    KPs_imageCoor_list=[] #[[stem],[leaf1],...,[ear]]
    stem_KP_list=[[]]
    leaves_KP_list=[]
    ear_KP_list=[[]]
    plant_KP_imageCoor_dict={} # KP: keyPoints
    plant_KP_virtualCoor_dict={} # KP: keyPoints
    
    line_i=0
    
    while line_i<dataIn_len:
        line=dataIn[line_i]
        if line.startswith('Scale'):
            line_i+=1
            line=dataIn[line_i]
            words=line.split(': ')
            try:
                scale=float(words[0].split(',')[1])
            except:
                raise Exception('Input key points file format error! Scale cannot be determined. Line: %d'%(line_i+1))
            if scale<=0:
                raise Exception('Scale setting error! scale<=0. Line: %d'%(line_i+1))
            line_i+=1
            line=dataIn[line_i]
        if line.startswith('<<'):
            if plant_info_str:
                stem_tip=KPs_imageCoor_list[-1][-1][0] # base of ear (or top-most leaf)
                stem_extendTip=KPs_imageCoor_list[0][-1][-1] # stem last point
                KPs_imageCoor_list[0][-1][-1]=stem_tip
                KPs_imageCoor_list[0][-1].append(stem_extendTip)
                plant_KP_imageCoor_dict[image_fn]=KPs_imageCoor_list
                stem_extendTip_KP=stem_KP_list[0][-1]
                if spikeInternodeLenZero and len(ear_KP_list[0])>0:
                    diff_vec=np.array(ear_KP_list[0][0]) - leaves_KP_list[-1][0]
                    ear_KP_list[0]=[ear_KP_list[0][i] - diff_vec for i in range(len(ear_KP_list[0]))]
                if ear_KP_list[0]:
                    stem_tip_KP=ear_KP_list[0][0]
                else:
                    stem_tip_KP=leaves_KP_list[-1][0]

                stem_KP_list[0][-1]=stem_tip_KP
                stem_KP_list[0].append(stem_extendTip_KP)
                if leafStraighten==1:
                    leaves_KP_list_new=[]
                    for leaf_i in leaves_KP_list:
                        vec_leaf_i=np.array(leaf_i[1])-leaf_i[0]
                        vec_leaf_i=vec_leaf_i/np.linalg.norm(vec_leaf_i)
                        len_leaf_i=[np.array(leaf_i[j])-leaf_i[j-1] for j in range(1,len(leaf_i))]
                        len_leaf_i=sum([np.linalg.norm(c) for c in len_leaf_i])
                        leaf_i_endP=leaf_i[0]+vec_leaf_i*len_leaf_i
                        leaves_KP_list_new.append([leaf_i[0],list(leaf_i_endP)])
                    leaves_KP_list=leaves_KP_list_new
                plant_KP_virtualCoor_dict[plant_info_str]={'stem':stem_KP_list,'leaf':leaves_KP_list,'ear':ear_KP_list}
                plant_info_str=''
            if keyWord in line.split(' ')[0]:
                image_fn=line[2:]
                plant_info_list=image_fn.split('.')[0].split(' ')
                cultivarID,plantID=plant_info_list[0].split('plant') # plant_info_list[0]=2-N-plant1,...
                cultivarID=cultivarID[0:-1]
                tillerID=plant_info_list[1][1:-1] # plant_info_list[1]=(1), (2), ...
                try:
                    int(tillerID)
                except:
                    print('Image file name: ',image_fn,'. Extracted tillerID: ',tillerID+'.')
                    raise Exception('ERROR: tillerID NOT a integer! Image file name error?')
                plant_info_str=cultivarID + '*' + plantID + '*' + tillerID
                KPs_imageCoor_list=[]
                stem_KP_list=[[]]
                leaves_KP_list=[]
                ear_KP_list=[[]]
                line_i += 1
                line=dataIn[line_i]
        try:
            int(line[0])
            data_line=True
        except:
            data_line=False
        if data_line and plant_info_str:
            words=line.split('\t')
            data_type=words[-1]
            #if data_type.startswith('Image_Crop_Coordinates'):
            #    [y1,y2,x1,x2]=re.split('[,;]+',words[0])
            if data_type.startswith('Stem'):
                stem_coor_list=[c.split(',') for c in words[0:-1]]
                stem_coor_list=[[float(c[0]),img_H-float(c[1])+1] for c in stem_coor_list]
                KPs_imageCoor_list.append([stem_coor_list])
                num_stem_KP=len(stem_coor_list)
                #z_direct=np.array(stem_coor_list[(num_stem_KP-1)//2+1]) - stem_coor_list[(num_stem_KP-1)//2]
                z_direct=np.array(stem_coor_list[-1]) - stem_coor_list[0]
                z_direct=z_direct/np.linalg.norm(z_direct)
                stem_coor_new=np.matrix(np.array(stem_coor_list)-stem_coor_list[0]).T
                stem_coor_new=[[z_direct[1],-z_direct[0]],[z_direct[0],z_direct[1]]]*stem_coor_new/scale
                
                stem_coor_new=stem_coor_new.T.tolist()
                stem_KP_list[-1]=stem_coor_new
            if data_type.startswith('Leaf'):
                leaves_KP_list.append([])
                leaf_coor_list=[c.split(',') for c in words[0:-1]]
                leaf_coor_list=[[float(c[0]),img_H-float(c[1])+1] for c in leaf_coor_list]
                KPs_imageCoor_list.append([leaf_coor_list])
                num_leaf_KP=len(leaf_coor_list)
                leaf_coor_new=np.matrix(np.array(leaf_coor_list)-stem_coor_list[0]).T
                leaf_coor_new=[[z_direct[1],-z_direct[0]],[z_direct[0],z_direct[1]]]*leaf_coor_new/scale
                leaf_coor_new=leaf_coor_new.T.tolist()
                leaves_KP_list[-1]=leaf_coor_new
            if data_type.startswith('Panicle'):
                ear_coor_list=[c.split(',') for c in words[0:-1]]
                ear_coor_list=[[float(c[0]),img_H-float(c[1])+1] for c in ear_coor_list]
                KPs_imageCoor_list.append([ear_coor_list])
                num_ear_KP=len(ear_coor_list)
                ear_coor_new=np.matrix(np.array(ear_coor_list)-stem_coor_list[0]).T
                ear_coor_new=[[z_direct[1],-z_direct[0]],[z_direct[0],z_direct[1]]]*ear_coor_new/scale
                ear_coor_new=ear_coor_new.T.tolist()
                ear_KP_list[-1]=ear_coor_new
        line_i+=1
    if plant_info_str:        
        stem_tip=KPs_imageCoor_list[-1][-1][0] # base of ear (or top-most leaf)
        stem_extendTip=KPs_imageCoor_list[0][-1][-1] # stem last point
        KPs_imageCoor_list[0][-1][-1]=stem_tip
        KPs_imageCoor_list[0][-1].append(stem_extendTip)
        plant_KP_imageCoor_dict[image_fn]=KPs_imageCoor_list
        stem_extendTip_KP=stem_KP_list[0][-1]
        if spikeInternodeLenZero and len(ear_KP_list[0])>0:
            diff_vec=np.array(ear_KP_list[0][0]) - leaves_KP_list[-1][0]
            ear_KP_list[0]=[ear_KP_list[0][i] - diff_vec for i in range(len(ear_KP_list[0]))]
        if ear_KP_list[0] and not spikeInternodeLenZero:
            stem_tip_KP=ear_KP_list[0][0]
        else:
            stem_tip_KP=leaves_KP_list[-1][0]
        stem_KP_list[0][-1]=stem_tip_KP
        stem_KP_list[0].append(stem_extendTip_KP)
        if leafStraighten==1:
            leaves_KP_list_new=[]
            for leaf_i in leaves_KP_list:
                vec_leaf_i=np.array(leaf_i[1])-leaf_i[0]
                vec_leaf_i=vec_leaf_i/np.linalg.norm(vec_leaf_i)
                len_leaf_i=[np.array(leaf_i[j])-leaf_i[j-1] for j in range(1,len(leaf_i))]
                len_leaf_i=sum([np.linalg.norm(c) for c in len_leaf_i])
                leaf_i_endP=leaf_i[0]+vec_leaf_i*len_leaf_i
                leaves_KP_list_new.append([leaf_i[0],list(leaf_i_endP)])
            leaves_KP_list=leaves_KP_list_new
        plant_KP_virtualCoor_dict[plant_info_str]={'stem':stem_KP_list,'leaf':leaves_KP_list,'ear':ear_KP_list}
    return plant_KP_imageCoor_dict,plant_KP_virtualCoor_dict

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

def buildLeaf(barcode,LNC_Kt_Kr,organ_PT_coords,stem_seg_base,stem_seg_tip,widLen_relation_coef,leafWidFun,rollRad,sheathRollRad,transLen,ptDense,halfSplit,randZRot_max,tortAng,tortAng_start_pos,tortAng_end_pos,leafAng_coef): # z axis stem:[stem_seg_base stem_seg_tip]
    ########### barcode=[plant_i,tiller_i,position_i,organCode]
    # para. consistency check
    tortAng_num=len(tortAng)
    if tortAng_num!=len(tortAng_start_pos) or tortAng_num!=len(tortAng_end_pos):
        raise Exception('Inconsistency of tortAng and their begin_pos and end_pos number!')
    for i in range(len(tortAng_start_pos)):
        if tortAng_start_pos[i]>tortAng_end_pos[i]:
            raise Exception('Inconsistency of tortAng begin_pos and end_pos!')
    # leaf segments and length
    leaf_segments_num=len(organ_PT_coords)-1
    leaf_segments_len=[np.linalg.norm(organ_PT_coords[i+1]-organ_PT_coords[i]) for i in range(leaf_segments_num)]
    leaf_len_max=sum(leaf_segments_len)
    # leafWidPattern
    leafWidMax=widLen_relation_coef[0]*leaf_len_max**2 + widLen_relation_coef[1]*leaf_len_max + widLen_relation_coef[2]
    # coord. transformation
    y_axis_new=np.array(stem_seg_tip) - stem_seg_base
    y_axis_new=y_axis_new/np.linalg.norm(y_axis_new)
    x_axis_new=[y_axis_new[2],y_axis_new[1],-y_axis_new[0]]
    z_axis_new=[0,-1,0]
    coord_transform_mat=np.array([x_axis_new,y_axis_new,z_axis_new])
    organ_PT_coords_new=np.dot(organ_PT_coords-stem_seg_base,coord_transform_mat.T)
    # adjusting leaf base angle
    if abs(leafAng_coef-1)>0.001:
        leaf_baseVec=organ_PT_coords_new[1] - organ_PT_coords_new[0]
        leaf_baseVec=leaf_baseVec/np.linalg.norm(leaf_baseVec)
        leaf_baseAng=math.acos(np.dot([0,1,0],leaf_baseVec))
        leaf_baseAng_new=leaf_baseAng*leafAng_coef
        leaf_baseAng_delta=leaf_baseAng_new-leaf_baseAng
        if leaf_baseVec[0]<0: # x < 0
            leaf_baseAng_delta=-leaf_baseAng_delta
        sin_delta=math.sin(leaf_baseAng_delta)
        cos_delta=math.cos(leaf_baseAng_delta)
        z_rot_mat=np.array([[cos_delta,-sin_delta,0],[sin_delta,cos_delta,0],[0,0,1]])
        organ_PT_coords_new=np.dot(organ_PT_coords_new-organ_PT_coords_new[0],z_rot_mat)+organ_PT_coords_new[0]
    organ_PT_coords_new_dense=[organ_PT_coords_new[0].tolist()] # leaf vein dense pt list
    t=[0] # cummulate leaf len list
    for i in range(leaf_segments_num):
        seg_i_vec=organ_PT_coords_new[i+1]-organ_PT_coords_new[i]
        seg_i_vec_len=np.linalg.norm(seg_i_vec)
        if seg_i_vec_len>0:
            seg_i_vec=seg_i_vec/seg_i_vec_len
        else:
            seg_i_vec=[0,0,0]
        seg_i_num_temp=leaf_segments_len[i]//ptDense
        for seg_i_j in range(1,int(seg_i_num_temp)):
            new_pt_temp=organ_PT_coords_new[i]+seg_i_vec*ptDense*seg_i_j
            organ_PT_coords_new_dense.append(new_pt_temp.tolist())
            t.append(t[-1]+ptDense)
        organ_PT_coords_new_dense.append(organ_PT_coords_new[i+1].tolist())
        t.append(sum(leaf_segments_len[0:i+1]))
    # leaf z axis rotation (y axis rotation in new coord.)
    randZRot=randZRot_max*2*(random.random()-0.5)
    rotation_matrix=np.array([[math.cos(randZRot),0,math.sin(randZRot)],[0,1,0],[-math.sin(randZRot),0,math.cos(randZRot)]])
    # 3D (curvature + roll) leaf in new coord.
    leaf_wid_list = leafWidMax*leafWidFun
    x2 = np.linspace(0,leaf_len_max,len(leafWidFun))
    y2 = leaf_wid_list/2
    interpolant = PchipInterpolator(x2,y2)
    p = interpolant(t) # leaf width list
    
    if rollRad==0:
        rollRad_y=100
    else:
        rollRad_y=rollRad
    
    tortAng_old=0
    tot_leafPts=len(t)
    row = np.array([0.0]*(2*halfSplit+1)*3*tot_leafPts).reshape((2*halfSplit+1,3,tot_leafPts))
    # pts setting: base -> tip
    transLen=min(max(0,transLen),leaf_len_max)
    for pt_i in range(tot_leafPts):
        if pt_i==tot_leafPts-1:
            x_axis_new2=np.array(organ_PT_coords_new_dense[pt_i])-organ_PT_coords_new_dense[pt_i-1]
            #print('***',organ_PT_coords_new_dense[pt_i],organ_PT_coords_new_dense[pt_i-1])
        else:
            x_axis_new2=np.array(organ_PT_coords_new_dense[pt_i+1])-organ_PT_coords_new_dense[pt_i]
            #print('***',organ_PT_coords_new_dense[pt_i+1],organ_PT_coords_new_dense[pt_i])
        x_axis_new2_len=np.linalg.norm(x_axis_new2)
        if x_axis_new2_len>0:
            x_axis_new2=x_axis_new2/np.linalg.norm(x_axis_new2)
        else:
            x_axis_new2=[1,0,0]
        y_axis_new2=[-x_axis_new2[1],x_axis_new2[0],0]
        z_axis_new2=[0,0,1]
        
        if pt_i < transLen//ptDense:
            R1=sheathRollRad*(transLen/ptDense-pt_i)/(transLen/ptDense)+rollRad_y*pt_i/(transLen/ptDense)
        else:
            R1=rollRad_y
        alpha=math.pi/2
        tortAng_now=tortAng_old
        for pos_i in range(len(tortAng_start_pos)):
            #print(pt_i,t[pt_i]/leaf_len_max,tortAng_start_pos[pos_i],tortAng_end_pos[pos_i],tortAng_now)
            if tortAng_start_pos[pos_i]<=t[pt_i]/leaf_len_max and tortAng_end_pos[pos_i]>=t[pt_i]/leaf_len_max:
                tortAng_now=sum(tortAng[0:pos_i])+tortAng[pos_i]*(t[pt_i]/leaf_len_max-tortAng_start_pos[pos_i])/(tortAng_end_pos[pos_i]-tortAng_start_pos[pos_i]) # leaf face direction alteration radian along leaf length
                break
        tortAng_old=tortAng_now
        # leaf twist matrix
        r_mat=np.array([[1,0,0],[0,math.cos(tortAng_now),-math.sin(tortAng_now)],[0,math.sin(tortAng_now),math.cos(tortAng_now)]])
        pt_coord_list=[]
        for j in range(2*halfSplit+1):
            alpha_new=alpha-p[pt_i]*(halfSplit-j)/halfSplit/R1
            #print(alpha,p[pt_i],halfSplit,R1)
            temp_coord=[0,R1]+R1*np.array([-math.cos(alpha_new),-math.sin(alpha_new)])
            pt_coord=[0,temp_coord[1],temp_coord[0]]
            pt_coord=np.dot(pt_coord,r_mat) # leaf twist
            pt_coord=np.dot(pt_coord,np.array([x_axis_new2,y_axis_new2,z_axis_new2])) + organ_PT_coords_new_dense[pt_i]
            row[j,:,pt_i]=pt_coord
    M_leaf=np.zeros([0,9])
    # triangularization: base -> tip
    for i in range(1,len(t)):
        pt_tri=np.zeros([4*halfSplit,9])
        for j in range(1,2*halfSplit+1):
            pt_tri[2*j-2,:]=np.hstack((row[j-1,0:3,i-1] , row[j-1,0:3,i] , row[j,0:3,i-1])) # tri_left_up
            pt_tri[2*j-1,:]=np.hstack((row[j-1,0:3,i] , row[j,0:3,i] , row[j,0:3,i-1])) # tri_right_up
        M_leaf=np.vstack((M_leaf,pt_tri))
        
    M_leaf_temp=M_leaf
    triangle_num=np.size(M_leaf,0)
    M_leaf=np.zeros([triangle_num,17])
    M_leaf[:,0:5]=np.tile(barcode,(triangle_num,1))
    M_leaf[:,5:8]=np.dot(M_leaf_temp[:,0:3],rotation_matrix)
    M_leaf[:,8:11]=np.dot(M_leaf_temp[:,3:6],rotation_matrix)
    M_leaf[:,11:14]=np.dot(M_leaf_temp[:,6:9],rotation_matrix)
    M_leaf[:,14:17]=np.tile(LNC_Kt_Kr,(triangle_num,1))
    # converting back to original old coords. (organ_PT_coords)
    off_set3=np.tile(stem_seg_base,(np.size(M_leaf,0),1)) 
    M_leaf[:,5:8]=off_set3 + np.dot(M_leaf[:,5:8],coord_transform_mat)
    M_leaf[:,8:11]=off_set3 + np.dot(M_leaf[:,8:11],coord_transform_mat)
    M_leaf[:,11:14]=off_set3 + np.dot(M_leaf[:,11:14],coord_transform_mat)
    return M_leaf

def buildSpikeletCover(barcode,LNC_Kt_Kr,awnLen,awnThick,grainLen,grainThick1,grainThick2,awnAng,allowedLayerNum):
    M_spikelet_cover=np.zeros([0,17])
    
    random_extent0=0.2
    random_extent1=0

    y_in_order=[1,1/2,-1/2,-1]
    x_in_order=[0,-math.sqrt(3)*1/2,-math.sqrt(3)*1/2,0]
    if allowedLayerNum==2:
        radi_foldChange=[1/2,1,0]
    elif allowedLayerNum==3:
        radi_foldChange=[1/2,1,3/5,0]
    elif allowedLayerNum==4:
        radi_foldChange=[1/2,4/5,1,1/2,0]
    else:
        raise Exception('Illegal allowed layer number for building spikelet cover!')
    for layer in range(allowedLayerNum): # 
        radi_x0=grainThick1/2*radi_foldChange[layer] # larger
        radi_y0=grainThick2/2*radi_foldChange[layer]
        radi_x1=grainThick1/2*radi_foldChange[layer+1] # larger
        radi_y1=grainThick2/2*radi_foldChange[layer+1]
        x0=0
        y0=0
        z0=grainLen*layer/allowedLayerNum
        x1=0
        y1=0
        z1=grainLen*(layer+1)/allowedLayerNum
        for i in range(3):
            if layer!=allowedLayerNum-1: # layer!=0 and allowedLayerNum-1
                p_s_01 = [radi_x0*x_in_order[i] +   x0, radi_y0*y_in_order[i] +   y0, z0]
                p_s_02 = [radi_x0*x_in_order[i+1] + x0, radi_y0*y_in_order[i+1] + y0, z0]
                p_s_11 = [radi_x1*x_in_order[i] +   x1, radi_y1*y_in_order[i] +   y1, z1]
                p_s_12 = [radi_x1*x_in_order[i+1] + x1, radi_y1*y_in_order[i+1] + y1, z1]
                M_patch = [barcode + p_s_01 + p_s_02 + p_s_11 + LNC_Kt_Kr[0:3]] + [barcode + p_s_02 + p_s_11 + p_s_12 + LNC_Kt_Kr[0:3]]
            else: # last layer
                p_s_01 = [radi_x0*x_in_order[i] +   x0, radi_y0*y_in_order[i] +   y0, z0]
                p_s_02 = [radi_x0*x_in_order[i+1] + x0, radi_y0*y_in_order[i+1] + y0, z0]
                p_s_11 = [x1, y1, z1]
                M_patch = [barcode + p_s_01 + p_s_02 + p_s_11 + LNC_Kt_Kr[0:3]]
            M_spikelet_cover = np.vstack((M_spikelet_cover,M_patch))
    # awn
    barcode_new=copy.deepcopy(barcode)
    barcode_new[2]=-5 # branch: -2; grain: -3; glume: -4; awn: -5
    if awnLen>0:
        x_in_order2=[0,-0.577,0.577,0] # 1/math.sqrt(3)
        y_in_order2=[2/3,-1/3,-1/3,2/3]
        radi_x=np.array(awnThick)/2
        radi_y=np.array(awnThick)/2
        x0=0-radi_x[0]
        y0=0
        z0=grainLen*0.95 # base of awn within grain
        awnAng_new=[0,0]
        awnAng_abs0=abs(awnAng[0])*10+0.0000001
        awnAng_abs1=abs(awnAng[1])*10+0.0000001
        awnAng_new[0]=-awnAng[0]#+math.sin(awnAng_abs0)/awnAng_abs0*(random.random()-0.5)*2*random_extent0 #*(1+(random.random()-0.5)*2*random_extent)
        awnAng_new[1]=-awnAng[1]#+math.sin(awnAng_abs1)/awnAng_abs1*(random.random()-0.5)*2*random_extent0 #*(1+(random.random()-0.5)*2*random_extent)
        # awnAng shrink to z axis in front view of a spikelet
        cos_yRotAng=math.cos(awnAng_new[0])
        sin_yRotAng=math.sin(awnAng_new[0])
        y_rot_mat=np.array([[cos_yRotAng,0,sin_yRotAng],[0,1,0],[-sin_yRotAng,0,cos_yRotAng]])
        # awnAng shrink to z axis in side view of a spikelet
        cos_xRotAng=math.cos(awnAng_new[1])
        sin_xRotAng=math.sin(awnAng_new[1])
        x_rot_mat=np.array([[1,0,0],[0,cos_xRotAng,-sin_xRotAng],[0,sin_xRotAng,cos_xRotAng]])
        # random small pertubation of awn direction around spikeletCover z axis
        zRotAng=math.pi*(random.random()-0.5)*2*random_extent1
        cos_zRotAng=math.cos(zRotAng)
        sin_zRotAng=math.sin(zRotAng)
        z_rot_mat=np.array([[cos_zRotAng,-sin_zRotAng,0],[sin_zRotAng,cos_zRotAng,0],[0,0,1]])
        x2=0
        y2=0
        z2=awnLen+grainLen*0.05
        [x2,y2,z2]=np.dot(np.dot(np.dot([x2,y2,z2],y_rot_mat),x_rot_mat),z_rot_mat)+[x0,y0,z0]
        [x1,y1,z1]=(np.array([x2,y2,z2])+[x0,y0,z0])/2
        for i in range(3):
            p_s_01 = [radi_x[0]*x_in_order2[i] +   x0, radi_y[0]*y_in_order2[i] +   y0, z0]
            p_s_02 = [radi_x[0]*x_in_order2[i+1] + x0, radi_y[0]*y_in_order2[i+1] + y0, z0]
            p_s_11 = [radi_x[1]*x_in_order2[i] +   x1, radi_y[1]*y_in_order2[i] +   y1, z1]
            p_s_12 = [radi_x[1]*x_in_order2[i+1] + x1, radi_y[1]*y_in_order2[i+1] + y1, z1]
            M_patch = [barcode_new + p_s_01 + p_s_02 + p_s_11 + LNC_Kt_Kr[3:]]+[barcode_new + p_s_02 + p_s_11 + p_s_12 + LNC_Kt_Kr[3:]]
            M_spikelet_cover = np.vstack((M_spikelet_cover,M_patch))
        for i in range(3):
            p_s_01 = [radi_x[1]*x_in_order2[i] +   x1, radi_y[1]*y_in_order2[i] +   y1, z1]
            p_s_02 = [radi_x[1]*x_in_order2[i+1] + x1, radi_y[1]*y_in_order2[i+1] + y1, z1]
            p_s_11 = [radi_x[2]*x_in_order2[i] +   x2, radi_y[2]*y_in_order2[i] +   y2, z2]
            p_s_12 = [radi_x[2]*x_in_order2[i+1] + x2, radi_y[2]*y_in_order2[i+1] + y2, z2]
            M_patch = [barcode_new + p_s_01 + p_s_02 + p_s_11 + LNC_Kt_Kr[3:]]+[barcode_new + p_s_02 + p_s_11 + p_s_12 + LNC_Kt_Kr[3:]]
            M_spikelet_cover = np.vstack((M_spikelet_cover,M_patch))
    return M_spikelet_cover
    
    
def buildSpikelet(barcode,LNC_Kt_Kr,spikeletBranchLen,spikeletBranchThick,awnLen_lemma,awnThick,grainLen,grainWid,grainThick,grainLen2,grainWid2,grainThick2,awnAng_coef,glumeAng,allowedLayerNum,glumeAwnLen,spikeletAng,LEMMA_ANG):
    # Height of the center spikelet (spikelet2)
    centerSpikeletH=0.2
    
    barcode_new=copy.deepcopy(barcode)
    M_spikelet_all=np.zeros([0,17])
    ### spikelet branch
    barcode_new[2]=-2 # branch: -2; grain: -3; glume: -4; awn: -5
    x0 = 0
    y0 = 0
    z0 = 0
    dem0 = spikeletBranchThick
    radi0 = dem0/2
    
    x1 = 0
    y1 = 0
    z1 = spikeletBranchLen
    dem1 = dem0*1
    radi1 = dem1/2
    
    x_in_order=[1,1/2,-1/2,-1,-1/2,1/2,1]
    y_in_order=[0,math.sqrt(3)*1/2,math.sqrt(3)*1/2,0,-math.sqrt(3)*1/2,-math.sqrt(3)*1/2,0]
    
    for i in range(6):
        p_s_01 = [radi0*x_in_order[i] +   x0, radi0*y_in_order[i] +   y0, z0]
        p_s_02 = [radi0*x_in_order[i+1] + x0, radi0*y_in_order[i+1] + y0, z0]
        p_s_11 = [radi1*x_in_order[i] +   x1, radi1*y_in_order[i] +   y1, z1]
        p_s_12 = [radi1*x_in_order[i+1] + x1, radi1*y_in_order[i+1] + y1, z1]
        M_spikeletBranch = [barcode_new + p_s_01 + p_s_02 + p_s_11 + LNC_Kt_Kr[0:3]] + [barcode_new + p_s_02 + p_s_11 + p_s_12 + LNC_Kt_Kr[0:3]]
        M_spikelet_all = np.vstack((M_spikelet_all,M_spikeletBranch))
    baseCoord_spikelets=[x1,y1,z1*1/3]
    # calculate lemma1 angle
#     b_len=centerSpikeletH+grainLen/2
#     c_len=grainLen/2
#     a_len=grainThick1*1 # estimated dist between center of spikelet1 and spikelet2
#     cos_a_lemma1=(b_len**2+c_len**2-a_len**2)/(2*b_len*c_len)
#     sin_a_lemma1=math.sqrt(1-cos_a_lemma1**2)
    #print(grainThick1,grainLen,centerSpikeletH)
    if LEMMA_ANG==-100:
        lemmaAng=1*math.atan(grainWid/2/math.sqrt(centerSpikeletH**2+grainLen*centerSpikeletH))
    else:
        lemmaAng=LEMMA_ANG
    sin_a_lemma1=math.sin(lemmaAng)
    cos_a_lemma1=math.cos(lemmaAng)
    if glumeAng<lemmaAng:
        print('WARNING: Glume angle < lemma angle, reset it as %f'%lemmaAng)
        glumeAng=lemmaAng
    cos_a_glume=math.cos(glumeAng)
    sin_a_glume=math.sin(glumeAng)
    ### glume
    barcode_new[2]=-4 # branch: -2; grain: -3; glume: -4; awn: -5
#     awnThick_glume=np.array(awnThick)/2
#     awnAng_glume=[0,0]
#     awnAng_glume[0]=glumeAng*awnAng_coef[0]
#     awnAng_glume[1]=spikeletAng*awnAng_coef[1]
#     M_glume1 = buildSpikeletCover(barcode_new,LNC_Kt_Kr,glumeAwnLen,awnThick_glume,grainLen,grainThick1,grainThick2,awnAng_glume,allowedLayerNum)
#     # glumeAng (y rotate)
#     triangle_num=np.size(M_glume1,0)
#     y_rotate1=np.array([[cos_a_glume,0,sin_a_glume],[0,1,0],[-sin_a_glume,0,cos_a_glume]])
    baseCoord_left=[-grainWid2/4,0,z1*1/3]
    baseCoord_right=[-baseCoord_left[0],0,z1*1/3]
#     off_set=np.tile(baseCoord_left,(triangle_num,1)) 
#     M_glume1[:,5:8]=off_set + np.dot(M_glume1[:,5:8],y_rotate1)
#     M_glume1[:,8:11]=off_set + np.dot(M_glume1[:,8:11],y_rotate1)
#     M_glume1[:,11:14]=off_set + np.dot(M_glume1[:,11:14],y_rotate1)
#     M_spikelet_all = np.vstack((M_spikelet_all,M_glume1))
#     M_glume2 = buildSpikeletCover(barcode_new,LNC_Kt_Kr,glumeAwnLen,awnThick_glume,grainLen,grainThick1,grainThick2,awnAng_glume,allowedLayerNum)
#     M_glume2[:,5:8]=off_set + np.dot(M_glume2[:,5:8],y_rotate1)
#     M_glume2[:,8:11]=off_set + np.dot(M_glume2[:,8:11],y_rotate1)
#     M_glume2[:,11:14]=off_set + np.dot(M_glume2[:,11:14],y_rotate1)
#     M_glume2[:,5]=-M_glume2[:,5]
#     M_glume2[:,8]=-M_glume2[:,8]
#     M_glume2[:,11]=-M_glume2[:,11]
#     M_spikelet_all = np.vstack((M_spikelet_all,M_glume2))
    ### lemma
    # lemma1 and lemma3
    awnThick_lemma=awnThick
    awnAng_lemma=[0,0]
    awnAng_lemma[0]=lemmaAng*awnAng_coef[0]
    awnAng_lemma[1]=spikeletAng*awnAng_coef[1]
    M_lemma1 = buildSpikeletCover(barcode_new,LNC_Kt_Kr,awnLen_lemma,awnThick_lemma,grainLen,grainWid,grainThick,awnAng_lemma,allowedLayerNum)
    y_rotate1=np.array([[cos_a_lemma1,0,sin_a_lemma1],[0,1,0],[-sin_a_lemma1,0,cos_a_lemma1]])
    triangle_num=np.size(M_lemma1,0)
    off_set=np.tile(baseCoord_left,(triangle_num,1)) 
    M_lemma1[:,5:8]=off_set + np.dot(M_lemma1[:,5:8],y_rotate1)
    M_lemma1[:,8:11]=off_set + np.dot(M_lemma1[:,8:11],y_rotate1)
    M_lemma1[:,11:14]=off_set + np.dot(M_lemma1[:,11:14],y_rotate1)
    M_spikelet_all = np.vstack((M_spikelet_all,M_lemma1))
    
    # M_lemma3=M_lemma
    M_lemma3 = buildSpikeletCover(barcode_new,LNC_Kt_Kr,awnLen_lemma,awnThick_lemma,grainLen,grainWid,grainThick,awnAng_lemma,allowedLayerNum)
    M_lemma3[:,5:8]=off_set + np.dot(M_lemma3[:,5:8],y_rotate1)
    M_lemma3[:,8:11]=off_set + np.dot(M_lemma3[:,8:11],y_rotate1)
    M_lemma3[:,11:14]=off_set + np.dot(M_lemma3[:,11:14],y_rotate1)
    M_lemma3[:,5]=-M_lemma3[:,5]
    M_lemma3[:,8]=-M_lemma3[:,8]
    M_lemma3[:,11]=-M_lemma3[:,11]
    M_spikelet_all = np.vstack((M_spikelet_all,M_lemma3))
    # lemma2
    baseCoord=np.array(baseCoord_spikelets) + [0,0,centerSpikeletH]
    sin_a_lemma2=sin_a_lemma1/20 # lemma of the center spikelet also rotate a little
    cos_a_lemma2=math.sqrt(1-sin_a_lemma2**2)
    lemma2Ang=math.asin(sin_a_lemma2)
    awnAng_lemma2=[0,0]
    awnAng_lemma2[0]=lemma2Ang*awnAng_coef[0] # lemma2Ang*(1+(random.random()-0.5)*2*0.05)
    awnAng_lemma2[1]=spikeletAng*awnAng_coef[1]
    M_lemma2 = buildSpikeletCover(barcode_new,LNC_Kt_Kr,awnLen_lemma,awnThick_lemma,grainLen2,grainWid2,grainThick2,awnAng_lemma2,allowedLayerNum)
    triangle_num=np.size(M_lemma2,0)
    off_set=np.tile(baseCoord,(triangle_num,1)) 
    y_rotate2=np.array([[cos_a_lemma2,0,sin_a_lemma2],[0,1,0],[-sin_a_lemma2,0,cos_a_lemma2]])
    M_lemma2[:,5:8]=off_set + np.dot(M_lemma2[:,5:8],y_rotate2)
    M_lemma2[:,8:11]=off_set + np.dot(M_lemma2[:,8:11],y_rotate2)
    M_lemma2[:,11:14]=off_set + np.dot(M_lemma2[:,11:14],y_rotate2)
    zRot=(random.random()>0.5)
    if zRot: # random mirror translation of lemma2
        M_lemma2[:,5]=-M_lemma2[:,5]
        M_lemma2[:,8]=-M_lemma2[:,8]
        M_lemma2[:,11]=-M_lemma2[:,11]
    M_spikelet_all = np.vstack((M_spikelet_all,M_lemma2))
    #print(glumeAng,awnAng_glume,awnAng_lemma,awnAng_lemma2)
    ### palea
    # palea1 and palea3
    awnLen_palea=0
    awnThick_palea=0
    awnAng_palea=[0,0]
    paleaAng=lemmaAng*0.9 # palea not completely coincides with lemma
    cos_a_palea1=math.cos(paleaAng)
    sin_a_palea1=math.sin(paleaAng)
    M_palea = buildSpikeletCover(barcode_new,LNC_Kt_Kr,awnLen_palea,awnThick_palea,grainLen,grainWid,grainThick,awnAng_palea,allowedLayerNum)
    triangle_num=np.size(M_palea,0)
    y_rotate3=np.array([[cos_a_palea1,0,-sin_a_palea1],[0,1,0],[sin_a_palea1,0,cos_a_palea1]])
    off_set=np.tile(baseCoord_right,(triangle_num,1)) 
    M_palea3 = copy.deepcopy(M_palea)
    M_palea3[:,5:8]=off_set + np.dot(M_palea[:,5:8],y_rotate3)
    M_palea3[:,8:11]=off_set + np.dot(M_palea[:,8:11],y_rotate3)
    M_palea3[:,11:14]=off_set + np.dot(M_palea[:,11:14],y_rotate3)
    M_spikelet_all = np.vstack((M_spikelet_all,M_palea3))
    M_palea1 = copy.deepcopy(M_palea3) # palea1 is mirror of palea3
    M_palea1[:,5]=-M_palea1[:,5]
    M_palea1[:,8]=-M_palea1[:,8]
    M_palea1[:,11]=-M_palea1[:,11]
    M_spikelet_all = np.vstack((M_spikelet_all,M_palea1))

    # palea2
    baseCoord=np.array(baseCoord_spikelets) + [0,0,centerSpikeletH]
    off_set=np.tile(baseCoord,(triangle_num,1)) 
    y_rotate2=np.array([[cos_a_lemma2,0,sin_a_lemma2],[0,1,0],[-sin_a_lemma2,0,cos_a_lemma2]])
    M_palea2 = buildSpikeletCover(barcode_new,LNC_Kt_Kr,awnLen_palea,awnThick_palea,grainLen2,grainWid2,grainThick2,awnAng_palea,allowedLayerNum)
    M_palea2[:,5:8]=off_set + np.dot(M_palea2[:,5:8],y_rotate2)
    M_palea2[:,8:11]=off_set + np.dot(M_palea2[:,8:11],y_rotate2)
    M_palea2[:,11:14]=off_set + np.dot(M_palea2[:,11:14],y_rotate2)
    if not zRot:
        M_palea2[:,5]=-M_palea2[:,5]
        M_palea2[:,8]=-M_palea2[:,8]
        M_palea2[:,11]=-M_palea2[:,11]
    M_spikelet_all = np.vstack((M_spikelet_all,M_palea2))
    # grain 2
#     M_grain21 = copy.deepcopy(M_palea)
#     M_grain21[:,2]=-3  # branch: -2; grain: -3; glume: -4; awn: -5
#     M_grain21[:,5:8]=off_set + M_palea[:,5:8]
#     M_grain21[:,8:11]=off_set + M_palea[:,8:11]
#     M_grain21[:,11:14]=off_set + M_palea[:,11:14]
#     M_spikelet_all = np.vstack((M_spikelet_all,M_grain21))
#     M_grain22 = copy.deepcopy(M_grain21)
#     M_grain22[:,5]=-M_grain22[:,5]
#     M_grain22[:,8]=-M_grain22[:,8]
#     M_grain22[:,11]=-M_grain22[:,11]
#     M_spikelet_all = np.vstack((M_spikelet_all,M_grain22))
    
    return M_spikelet_all
    
    
def buildEar(barcode,LNC_Kt_Kr,organ_PT_coords,spikeletNum,spikeletBranchLen,spikeletBranchThick,awnLen_list,awnThick,grainLen,grainWid,grainThick,grainLen2,grainWid2,grainThick2,grainShrinkCoef,awnAng_coef_bot,awnAng_coef_mid,awnAng_coef_top,glumeAng,lemmaAng,SPIKELET_ANG,zRotAng,earRachis_diameter,allowedLayerNum,glumeAwnLen):
    ########### barcode=[plant_i,tiller_i,position_i,organCode]
    # para. consistency check
    
    ### 3D reconstruction of spike rachis
    barcode_new=copy.deepcopy(barcode)
    organ_PT_coords_new=copy.deepcopy(organ_PT_coords)
    M_ear_all=np.zeros([0,17])
    ear_segments_num=len(organ_PT_coords_new)-1
    barcode_new[2]=-2 # branch: -2; grain: -3; glume: -4; awn: -5
    for pt_i in range(1,ear_segments_num):
        organ_base=organ_PT_coords_new[pt_i-1]
        organ_tip=organ_PT_coords_new[pt_i]
        M_ear = buildStem(barcode_new,LNC_Kt_Kr[0:3],organ_base,organ_tip,[earRachis_diameter,earRachis_diameter]) # plant, tiller, Position, organCode
        M_ear_all=np.vstack((M_ear_all,M_ear))
    pt_i=ear_segments_num
    organ_base=organ_PT_coords_new[pt_i-1]
    organ_vec=organ_PT_coords_new[pt_i]-organ_base
    organ_len=np.linalg.norm(organ_vec)
    organ_direction=(organ_vec)/organ_len
    organ_PT_coords_new[pt_i]=organ_base+organ_direction*(organ_len-grainShrinkCoef*grainLen)
    organ_tip=organ_PT_coords_new[pt_i]
    M_ear = buildStem(barcode_new,LNC_Kt_Kr[0:3],organ_base,organ_tip,[earRachis_diameter,earRachis_diameter*0.3]) # plant, tiller, Position, organCode
    M_ear_all=np.vstack((M_ear_all,M_ear))
    # ear segments and length
    ear_segments_len=[np.linalg.norm(organ_PT_coords_new[i+1]-organ_PT_coords_new[i]) for i in range(ear_segments_num)]
    ear_len_max=sum(ear_segments_len)
    spikeletDist=ear_len_max/spikeletNum
    ### 3D reconstruction of spikelets
    z_rotate_other=np.array([[math.cos(zRotAng),-math.sin(zRotAng),0],[math.sin(zRotAng),math.cos(zRotAng),0],[0,0,1]])
    zRotAng=zRotAng+math.pi/2
    z_rotate_tip=np.array([[math.cos(zRotAng),-math.sin(zRotAng),0],[math.sin(zRotAng),math.cos(zRotAng),0],[0,0,1]])
    z_rotate=z_rotate_other
    for sp_i in range(1,spikeletNum+1):
        grainLen_i=grainLen
        grainWid_i=grainWid
        grainThick_i=grainThick
        grainLen2_i=grainLen2
        grainWid2_i=grainWid2
        grainThick2_i=grainThick2
        sp_i_rank=sp_i/spikeletNum
        awnAng_coef_i=[1,1]
        if sp_i_rank<1/3:
            awnLen=awnLen_list[0]
            awnAng_coef_i[0]=awnAng_coef_bot[0]*(1-sp_i_rank*3)+awnAng_coef_mid[0]*sp_i_rank*3
            awnAng_coef_i[1]=awnAng_coef_bot[1]*(1-sp_i_rank*3)+awnAng_coef_mid[1]*sp_i_rank*3
        elif sp_i_rank<2/3:
            awnLen=awnLen_list[1]
            awnAng_coef_i=awnAng_coef_mid
        else:
            awnLen=awnLen_list[2]
            awnAng_coef_i[0]=awnAng_coef_mid[0]*(3-sp_i_rank*3)+awnAng_coef_top[0]*(sp_i_rank*3-2)
            awnAng_coef_i[1]=awnAng_coef_mid[1]*(3-sp_i_rank*3)+awnAng_coef_top[1]*(sp_i_rank*3-2)
            dist1=(spikeletNum-sp_i)*3/spikeletNum
            dist2=(sp_i-2*spikeletNum/3)*3/spikeletNum 
            grainLen_i=grainLen*dist1+grainShrinkCoef*grainLen*dist2
            grainWid_i=grainWid*dist1+grainShrinkCoef*grainWid*dist2
            grainThick_i=grainThick*dist1+grainShrinkCoef*grainThick*dist2
            grainLen2_i=grainLen2*dist1+grainShrinkCoef*grainLen2*dist2
            grainWid2_i=grainWid2*dist1+grainShrinkCoef*grainWid2*dist2
            grainThick2_i=grainThick2*dist1+grainShrinkCoef*grainThick2*dist2
        if SPIKELET_ANG==-100:
            spikeletAng=(-1)**sp_i*0.8*math.atan(grainThick_i/2/math.sqrt(spikeletBranchLen**2+grainLen_i*spikeletBranchLen)) # angle of spikelets to spike rachis
        else:
            spikeletAng=(-1)**sp_i*SPIKELET_ANG
        if sp_i==spikeletNum:
            spikeletAng=0  # the last spikelet is along the rachis direction
            z_rotate=z_rotate_tip
        dist_i=sp_i*spikeletDist
        for seg_i in range(ear_segments_num):
            seg_len_tot=sum(ear_segments_len[0:seg_i+1])
            if dist_i<=seg_len_tot+0.00001:
                z_vec=organ_PT_coords_new[seg_i+1]-organ_PT_coords_new[seg_i]
                z_axis_new=z_vec/np.linalg.norm(z_vec)
                Origin_coord=organ_PT_coords_new[seg_i]+(dist_i-sum(ear_segments_len[0:seg_i]))*z_axis_new
                break
        y_axis_new=[0,1,0]
        x_axis_new=[z_axis_new[2],z_axis_new[1],-z_axis_new[0]]
        coord_transform_mat=np.array([x_axis_new,y_axis_new,z_axis_new])
        
        if sp_i==spikeletNum:
            awnAng_coef_i=[1,1]
        #awnAng_coef_i[0]=min(1,awnAng_coef[0]*(2-sp_i_rank))
        #awnAng_coef_i[1]=min(1,awnAng_coef[1]*(2-sp_i_rank))
        M_spikelet = buildSpikelet(barcode,LNC_Kt_Kr[3:],spikeletBranchLen,spikeletBranchThick,awnLen,awnThick,grainLen_i,grainWid_i,grainThick_i,grainLen2_i,grainWid2_i,grainThick2_i,awnAng_coef_i,glumeAng,allowedLayerNum,glumeAwnLen,spikeletAng,lemmaAng)
        # spikeletAng (x rotate)
        x_rotate=np.array([[1,0,0],[0,math.cos(spikeletAng),-math.sin(spikeletAng)],[0,math.sin(spikeletAng),math.cos(spikeletAng)]])
        M_spikelet_temp=M_spikelet
        M_spikelet[:,5:8]=np.dot(M_spikelet_temp[:,5:8],x_rotate)
        M_spikelet[:,8:11]=np.dot(M_spikelet_temp[:,8:11],x_rotate)
        M_spikelet[:,11:14]=np.dot(M_spikelet_temp[:,11:14],x_rotate)
        # zRotAng (z rotate)
        M_spikelet[:,5:8]=np.dot(M_spikelet[:,5:8],z_rotate)
        M_spikelet[:,8:11]=np.dot(M_spikelet[:,8:11],z_rotate)
        M_spikelet[:,11:14]=np.dot(M_spikelet[:,11:14],z_rotate)
        # converting back to original old coords. (organ_PT_coords_new)
        triangle_num=np.size(M_spikelet,0)
        off_set=np.tile(Origin_coord,(triangle_num,1)) 
        M_spikelet[:,5:8]=off_set + np.dot(M_spikelet[:,5:8],coord_transform_mat)
        M_spikelet[:,8:11]=off_set + np.dot(M_spikelet[:,8:11],coord_transform_mat)
        M_spikelet[:,11:14]=off_set + np.dot(M_spikelet[:,11:14],coord_transform_mat)
        M_ear_all=np.vstack((M_ear_all,M_spikelet))
    return M_ear_all


def buildCanopy(M_multi_tillers,row_dist,sowing_width,plant_density,rowNum_simulate,rowLen_simulate,tillerAng,Kr_soil):
    
    TillerSize_list=list(map(lambda x:np.size(x,0),M_multi_tillers))
    maxTillerSize=max(TillerSize_list)
    
    plantsNum_oneRow=round(plant_density*rowLen_simulate)
    print(plantsNum_oneRow,maxTillerSize)
    plantNum_allRow=plantsNum_oneRow*rowNum_simulate
    M_canopy=np.zeros([plantNum_allRow*maxTillerSize+round(row_dist*(rowNum_simulate-1)*2*rowLen_simulate),17])
    print(np.size(M_canopy),np.size(M_canopy,1),np.size(M_canopy,0))
    tillerNum=len(M_multi_tillers)
    print('  In all ' + str(tillerNum) + ' tillers collected, building canopy with %d plants (rowLen:%dcm, rowNum:%d):'%(plantNum_allRow,rowLen_simulate,rowNum_simulate))
    start_pos=0
    for row_i in range(rowNum_simulate):
        print('    Row '+str(row_i+1)+'/'+str(rowNum_simulate)+'...')
        x_shiftDist_list=[]
        y_shiftDist_list=[]
        for plant_i in range(plantsNum_oneRow):
            x_shiftDist=row_i*row_dist+5*(random.random()-0.5)
            y_shiftDist=random.random()*rowLen_simulate
            zRotAng=math.pi*2*random.random() # tiller self rotate
            cos_zRotAng=math.cos(zRotAng)
            sin_zRotAng=math.sin(zRotAng)
            z_rot_mat=np.array([[cos_zRotAng,-sin_zRotAng,0],[sin_zRotAng,cos_zRotAng,0],[0,0,1]])
            tillerAng_new=tillerAng*random.random()
            cos_yRotAng=math.cos(tillerAng_new)
            sin_yRotAng=math.sin(tillerAng_new)
            y_rot_mat=np.array([[cos_yRotAng,0,sin_yRotAng],[0,1,0],[-sin_yRotAng,0,cos_yRotAng]])
            r_mat=np.dot(y_rot_mat,z_rot_mat)
            tiller_i=random.choice(range(tillerNum))
            xy_move=[x_shiftDist,y_shiftDist,0]
            triangle_num=np.size(M_multi_tillers[tiller_i],0)
            off_set=np.tile(xy_move,(triangle_num,1)) 
            M_tiller=np.zeros([triangle_num,17])
            M_tiller[:,0:5]=M_multi_tillers[tiller_i][:,0:5]
            M_tiller[:,14:]=M_multi_tillers[tiller_i][:,14:]
            M_tiller[:,5:8]=off_set + np.dot(M_multi_tillers[tiller_i][:,5:8],r_mat)
            M_tiller[:,8:11]=off_set + np.dot(M_multi_tillers[tiller_i][:,8:11],r_mat)
            M_tiller[:,11:14]=off_set + np.dot(M_multi_tillers[tiller_i][:,11:14],r_mat)
            TillerSize_current=TillerSize_list[tiller_i]
            current_pos=start_pos+TillerSize_current
            M_canopy[start_pos:current_pos,:]=M_tiller
            start_pos=current_pos
    xmin=0
    ymin=0
    zmin=0.01
    xmax=xmin+row_dist*(rowNum_simulate-1)
    ymax=ymin+rowLen_simulate
    # construct the soil horizontal plane
    for x_i in np.linspace(xmin,xmax,round(xmax-xmin)+1)[0:-1]:
        for y_i in np.linspace(ymin,ymax,round(ymax-ymin)+1)[0:-1]:
            p1=[x_i,y_i,zmin]
            p2=[x_i,y_i+1,zmin]
            p3=[x_i+1,y_i,zmin]
            p4=[x_i+1,y_i+1,zmin]
            M_soil = np.array([[0,0,-6,1,0] + p1 + p2 + p4 + [0,0,Kr_soil]] + [[0,0,-6,1,0] + p1 + p3 + p4 + [0,0,Kr_soil]])
            M_canopy[current_pos:current_pos+2,:]=M_soil
            current_pos+=2
    return M_canopy[0:current_pos,:]

def SPAD2Chl(SPAD):
    # gN/m2. -> ug/cm2. this equation combined two high-cited papers, the logic is SPAD -> LN (paper1); SPAD -> Chl (paper2); then LN -> Chl (solve eqn.)
    try:
        Chl=0.0599*math.exp(0.0493*SPAD)
    except:
        Chl=list(map(lambda x:0.0599*math.exp(0.0493*x),SPAD))
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

def TillerSkeletonPlot(plant_KP_virtualCoor_dict,figOutPrefix,showOn):
    xlow=-30
    xup=30
    ylow=-5
    yup=90
    
    plt.figure(figsize=(6*(xup-xlow)/(yup-ylow),6))
    ax=plt.axes()
    stem_PTs=plant_KP_virtualCoor_dict['stem'][0]
    stem_PTs_num=len(stem_PTs)
    for pt_i in range(1,stem_PTs_num-1): # stem_PTs_num-1: do not use the last point (extendTip)
        plt.plot([stem_PTs[pt_i-1][0],stem_PTs[pt_i][0]],[stem_PTs[pt_i-1][1],stem_PTs[pt_i][1]],color='k',marker='*',linestyle='-', linewidth=5, markersize=6)
    num_leaf=len(plant_KP_virtualCoor_dict['leaf'])
    for leaf_i in range(num_leaf):
        leaf_PTs=plant_KP_virtualCoor_dict['leaf'][leaf_i]
        leaf_PTs_num=len(leaf_PTs)
        for pt_i in range(1,leaf_PTs_num):
            plt.plot([leaf_PTs[pt_i-1][0],leaf_PTs[pt_i][0]],[leaf_PTs[pt_i-1][1],leaf_PTs[pt_i][1]],color='g',marker='o',linestyle='-', linewidth=3, markersize=6)
    ear_PTs=plant_KP_virtualCoor_dict['ear'][0]
    ear_PTs_num=len(ear_PTs)
    for pt_i in range(1,ear_PTs_num): # stem_PTs_num-1: do not use the last point (extendTip)
        plt.plot([ear_PTs[pt_i-1][0],ear_PTs[pt_i][0]],[ear_PTs[pt_i-1][1],ear_PTs[pt_i][1]],color='y',marker='s',linestyle='-', linewidth=8, markersize=6)
    plt.axis([xlow,xup,ylow,yup])
    #plt.axis('equal')
    plt.gca().set_aspect('equal', adjustable='box')
    if showOn:
        plt.show()
    plt.savefig(figOutPrefix+'.png',format="png",dpi=600)
    plt.close()
    
def awnAngCoef_Fitter(glumeAng,spikeletAng,awnAng1,awnAng2):
    awnAngCoef1=0
    awnAngCoef2=0
    fit_error=fit_residuals([awnAngCoef1,awnAngCoef2],glumeAng,spikeletAng,awnAng1,awnAng2)
    itera_test=0
    while itera_test < 50:
        #print('            Awn angle fitting: iteration %d, fit_error %.3f'%(itera_test,fit_error))
        awnAngCoef1_ini=random.uniform(-1,1)
        awnAngCoef2_ini=random.uniform(-1,1)
        paras0=np.array([awnAngCoef1_ini,awnAngCoef2_ini]) # theta,Amax,alpha,Rd
        popt= fmin_slsqp(fit_residuals, paras0, bounds=[(-1,1),(-1,1)], args=(glumeAng,spikeletAng,awnAng1,awnAng2),iprint = 0)
        popt=popt[0:]
        fit_error_new=fit_residuals(popt[0:],glumeAng,spikeletAng,awnAng1,awnAng2)
        if fit_error_new<fit_error:
            fit_error=fit_error_new
            awnAngCoef1,awnAngCoef2=popt
        itera_test=itera_test+1
    print('            Awn angle fitting done! awnAngCoef: %.2f, %.2f. Fit_error: %.3f'%(awnAngCoef1,awnAngCoef2,fit_error))
    return popt

def fit_residuals(paras,glumeAng,spikeletAng,awnAng1,awnAng2):
    awnAngCoef1,awnAngCoef2=paras
    z0=[0,0,1]
    x0=[1,0,0]
    ###spikelet-rachis angle
    cos_1=math.cos(-spikeletAng)
    sin_1=math.sin(-spikeletAng)
    #x_rot_mat=np.array([[1,0,0],[0,cos_1,-sin_1],[0,sin_1,cos_1]])
    z1=[0,sin_1,cos_1]#np.dot(z0,x_rot_mat)
    ###glume-spikelet angle
    rotation_vec=np.cross(z1,x0)
    rotation_vec=rotation_vec/np.linalg.norm(rotation_vec)
    cos_a=math.cos(glumeAng*(1-awnAngCoef1))
    sin_a=math.sin(glumeAng*(1-awnAngCoef1))
    a=rotation_vec[0]
    b=rotation_vec[1]
    c=rotation_vec[2]
    row_1=[a**2+(1-a**2)*cos_a,a*b*(1-cos_a)+c*sin_a,a*c*(1-cos_a)-b*sin_a]
    row_2=[a*b*(1-cos_a)-c*sin_a,b**2+(1-b**2)*cos_a,b*c*(1-cos_a)+a*sin_a]
    row_3=[a*c*(1-cos_a)+b*sin_a,c*b*(1-cos_a)-a*sin_a,c**2+(1-c**2)*cos_a]
    r_mat=np.array([row_1,row_2,row_3])
    z2=np.dot(z1,r_mat)
    cos_2=math.cos(spikeletAng*awnAngCoef2)
    sin_2=math.sin(spikeletAng*awnAngCoef2)
    x_rot_mat=np.array([[1,0,0],[0,cos_2,-sin_2],[0,sin_2,cos_2]])
    z3=np.dot(z2,x_rot_mat)
    #print('z3:',z3)
    awnAng1_fit=math.atan(z3[0]/z3[2])
    awnAng2_fit=math.atan(-z3[1]/z3[2])
    fit_err=(awnAng1_fit-awnAng1)**2+(awnAng2_fit-awnAng2)**2
    return fit_err

if __name__=="__main__":
    start=time.time()
    
    op=createHelp()
    log2file=op.log
    if log2file:
        sys.stdout = Logger('./log_wheatCanopyConstructor2')
    width_format=94
    welcome_info='Welcome to use wheat canopy constructor'
    time_info=time.asctime( time.localtime(time.time()) )
    print(welcome_info.center(width_format, '-'))
    print(time_info.center(width_format, ' '))
    runID=op.runID
    random.seed(runID)
    img_H=op.img_H
    cultivarID=op.lineNA
    fnOut_str=op.file_out
    overWrite=op.ow
    fnOut_list=fnOut_str.split('\\')
    resultDir='\\'.join(fnOut_list[0:-1])
    tillerSkeletonDir=os.path.join(resultDir,'temp_tillerSkeleton')
    if not os.path.exists(resultDir):
        print('  Result directory not exists, creating it ...')
        try:
            os.mkdir(resultDir)
        except:
            raise Exception('Cannot create result directory: %s'%resultDir)
    if 'Run' not in fnOut_list[-1]:
        fnOut_list[-1]='Run'+str(runID)+'_'+fnOut_list[-1]
    fnOut='\\'.join(fnOut_list)
    fnOut_plant=fnOut[0:-4]+'_plant.txt'
    op.file_out=fnOut
    inputDir=op.dir_input
    withSpike=op.withSpike
    config_morphology_fn,config_physiology_fn=op.file_in2
    config_morphology_fn=os.path.join(inputDir,config_morphology_fn)
    config_physiology_fn=os.path.join(inputDir,config_physiology_fn)

    data = json.dumps(vars(op), indent=4,ensure_ascii=False, sort_keys=False,separators=(',', ':'))
    print('Running parameters:')
    print(data)
    
    if (not overWrite) and os.path.exists(fnOut):
        print('  Existing 3D canopy structure file "%s" found! Nothing more needs to do, exit directly.'%fnOut)
        sys.exit()
    
    earLen_coef=1
    spikeInternodeLenZero=0
    leafStraighten=0
    transLen=0
    sheathRollRad=0
    rollRad=0
    allowedLayerNum=2
    ptDense=1
    halfSplit=2
    
    # read in parameters from configure files
    for line in open(config_morphology_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='plant_density':
            plant_density=float(words[1])
        elif words[0]=='tillerAng':
            tillerAng=float(words[1])
        elif words[0]=='rowNum':
            rowNum_simulate=int(words[1])# rows, rows simulated for a canopy
        elif words[0]=='rowLen':
            rowLen_simulate=float(words[1]) # cm, length of one row simulated for a canopy
        elif words[0]=='row_dist':
            row_dist=float(words[1]) #20 # cm, distance between two growing rows
        elif words[0]=='sowing_width':
            sowing_width=float(words[1])#10 # cm, the width of wheat seeds spread in one row
        elif words[0]=='plantH_coef':
            plantH_coef=float(words[1])#plant height scaling
        elif words[0]=='earLen_coef':
            earLen_coef=float(words[1])#ear length scaling
        elif words[0]=='spikeInternodeLenZero':
            spikeInternodeLenZero=float(words[1])#set spike Internode Length to 0
        elif words[0]=='leafStraighten':
            leafStraighten=float(words[1])# straighten all leaves
        elif words[0]=='widLen_relation_coef':
            widLen_relation_coef=float(words[1])#plant height scaling
        elif words[0]=='widLen_relation_flag_coef':
            widLen_relation_flag_coef=np.array([float(c) for c in words[1:]])
        elif words[0]=='widFunc_flag':
            widFunc_flag=[float(c) for c in words[1:]]
        elif words[0]=='widLen_relation_other_coef':
            widLen_relation_other_coef=np.array([float(c) for c in words[1:]])
        elif words[0]=='widFunc_other':
            widFunc_other=[float(c) for c in words[1:]]
        elif words[0]=='widLen_relation_developing_coef':
            widLen_relation_developing_coef=np.array([float(c) for c in words[1:]])
        elif words[0]=='widFunc_developing':
            widFunc_developing=[float(c) for c in words[1:]]
        elif words[0]=='leafAng_coef':
            leafAng_coef=float(words[1])#leaf angle scaling
        elif words[0]=='allowedLayerNum':
            allowedLayerNum=int(words[1]) # Detailed layers for building glume, lemma and grain. Can only be 2,3 or 4.
        elif words[0]=='rollRad':
            rollRad=float(words[1]) # leaf rolling radius
        elif words[0]=='sheathRollRad':
            sheathRollRad=float(words[1]) # stem radius, useless
        elif words[0]=='transLen':
            transLen=float(words[1]) # useless
        elif words[0]=='ptDense':
            ptDense=float(words[1]) # cm, cut leaf length direction into pieces with a size of ptDense
        elif words[0]=='halfSplit':
            halfSplit=int(words[1]) # cut leaf width direction into halfSplit*2 pieces
        elif words[0]=='randZRot_max':
            randZRot_max=float(words[1]) # random leaf rotating angle around z axis
        elif words[0]=='tortAng':
            tortAng=[float(c) for c in words[1:]] # leaf twisting angle
        elif words[0]=='tortAng_start_pos':
            tortAng_start_pos=[float(c) for c in words[1:]] # leaf twisting start position (relative to total leaf length)
        elif words[0]=='tortAng_end_pos':
            tortAng_end_pos=[float(c) for c in words[1:]] # leaf twisting end position (relative to total leaf length)
            
        elif words[0]=='stem_diameter':
            stem_diameter=np.array([float(c) for c in words[1:]])*1.047 # *1.047 to ensure Perimeter of hexagon equals to circle (area ratio 91%)
        elif words[0]=='stem_midH':
            stem_midH=float(words[1]) # cm, >30 belongs to upper, <30 belongs to lower
        elif words[0]=='spikeletDens':
            spikeletDense=float(words[1])
        elif words[0]=='grainShrinkCoef':
            grainShrinkCoef=float(words[1])
        elif words[0]=='earRachis_diameter':
            earRachis_diameter=float(words[1])
        elif words[0]=='spikeletBranchLen':
            spikeletBranchLen=float(words[1])
        elif words[0]=='spikeletBranchThick':
            spikeletBranchThick=float(words[1])
        elif words[0]=='allowedLayerNum':
            allowedLayerNum=int(words[1])
        elif words[0]=='grainEnlarge_coef':
            enlarge_coef=float(words[1])
        elif words[0]=='Side-grain length':
            grainLen=float(words[1])*enlarge_coef
        elif words[0]=='Side-grain width':
            grainWid=float(words[1])*enlarge_coef
        elif words[0]=='Mid-grain length':
            grainLen2=float(words[1])*enlarge_coef
        elif words[0]=='Mid-grain width':
            grainWid2=float(words[1])*enlarge_coef
        elif words[0]=='Tri-spikelet to rachis angle':
            SPIKELET_ANG=float(words[1])*math.pi/180
        elif words[0]=='Side-grain to mid-grain angle':
            LEMMA_ANG=float(words[1])*math.pi/180
            LEMMA_ANG=math.atan(math.tan(LEMMA_ANG)*math.cos(SPIKELET_ANG)) # convert project angle to plain angle 
        elif words[0]=='Glume to mid-grain angle':
            glumeAng=float(words[1])*math.pi/180
            glumeAng=math.atan(math.tan(glumeAng)*math.cos(SPIKELET_ANG)) # convert project angle to plain angle 
        elif words[0]=='Awn to rachis angle (front view)':
            awnAng1=float(words[1])*math.pi/180
        elif words[0]=='Awn to rachis angle (side view)':
            awnAng2=float(words[1])*math.pi/180
        elif words[0]=='Upper spike awn length':
            awnLen_upper=float(words[1])
        elif words[0]=='Middle spike awn length':
            awnLen_mid=float(words[1])
        elif words[0]=='Lower spike awn length':
            awnLen_lower=float(words[1])
        elif words[0]=='Awn base thickness':
            awnThick_base=float(words[1])
        elif words[0]=='Awn center thickness':
            awnThick_mid=float(words[1])
        elif words[0]=='Awn tip thickness':
            awnThick_tip=float(words[1])
    widLen_relation_flag_coef=widLen_relation_flag_coef*widLen_relation_coef
    widLen_relation_other_coef=widLen_relation_other_coef*widLen_relation_coef
    widLen_relation_developing_coef=widLen_relation_developing_coef*widLen_relation_coef
    for line in open(config_physiology_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='leaf_NC':
            leaf_LNC_distribution=np.array([float(c) for c in words[1:]])
            leaf_SPAD_distribution=9.70*leaf_LNC_distribution + 7.13 # LNC: %
        elif words[0]=='grain_NC':
            spikelet_LNC=float(words[1])
            spikelet_SPAD=9.70*spikelet_LNC + 7.13
        elif words[0]=='awn_NC':
            awn_LNC=float(words[1])
            awn_SPAD=9.70*awn_LNC + 7.13
        elif words[0]=='culm_NC':
            stem_LNC_distribution=np.array([float(c) for c in words[1:]]) # upper, lower
            stem_SPAD_distribution=9.70*stem_LNC_distribution + 7.13
        elif words[0]=='Kr_soil':
            Kr_soil=float(words[1]) # soil reflectance
    grainThick=grainWid
    grainThick2=grainWid2
    awnLen=[awnLen_lower, awnLen_mid, awnLen_upper]
    awnThick=[awnThick_base, awnThick_mid, awnThick_tip]
    glumeAwnLen=0
    
    spikeletAng=SPIKELET_ANG
    centerSpikeletH=0.2
    lemmaAng=LEMMA_ANG
    awnAng_coef_mid=awnAngCoef_Fitter(lemmaAng,spikeletAng,awnAng1,awnAng2) # ,awnAng1,awnAng2
    awnAng_coef_bot=awnAngCoef_Fitter(lemmaAng,spikeletAng,awnAng1/2,awnAng2*2) # awnAng1/2,awnAng2*2
    grainThick=(1+grainShrinkCoef)/2*grainThick
    grainLen=(1+grainShrinkCoef)/2*grainLen
    awnAng_coef_top=awnAngCoef_Fitter(lemmaAng,spikeletAng,awnAng1/4,awnAng2/4) # awnAng1/4,awnAng2/4
    
    
    zRotAng=90*math.pi/180 # rotation of spikelets around rachis, useless
    leaf_Chl_distribution=SPAD2Chl(leaf_SPAD_distribution)
    Kr_leaf_distribution=Chl2Kr(leaf_Chl_distribution)
    Kt_leaf_distribution=Chl2Kt(leaf_Chl_distribution)
    leaf_SPAD_mean=np.mean(leaf_SPAD_distribution[0:5])
    leaf_LNC_mean=np.mean(leaf_LNC_distribution[0:5])
    stem_Chl_distribution=SPAD2Chl(stem_SPAD_distribution)
    Kr_stem_distribution=Chl2Kr(stem_Chl_distribution)
    Kt_stem_distribution=[0,0]
    spikelet_Chl=SPAD2Chl(spikelet_SPAD)
    
    Kr_spikelet=Chl2Kr(spikelet_Chl)
    Kt_spikelet=Chl2Kt(spikelet_Chl)
    rachis_LNC=spikelet_LNC
    Kr_rachis=Kr_spikelet
    Kt_rachis=Kt_spikelet
    awn_Chl=SPAD2Chl(awn_SPAD)
    Kr_awn=Chl2Kr(awn_Chl)
    Kt_awn=Chl2Kt(awn_Chl)
    
    # read in structure of tillers
    dataIn=open(os.path.join(inputDir,op.file_in),'r').readlines()
    plant_KP_imageCoor_dict,plant_KP_virtualCoor_dict=inputFileProcessing(dataIn,img_H,cultivarID,spikeInternodeLenZero,leafStraighten)
    
    Plants_Tillers_info=list(plant_KP_virtualCoor_dict.keys())
    Plants_Tillers_info.sort(key=lambda x:int(x.split('*')[1])*100+int(x.split('*')[2]))
    M_multi_tillers=[]
    print('  Tillers in reconstructing...')
    for block_i in Plants_Tillers_info:
        if 0:
            if not os.path.exists(tillerSkeletonDir):
                print('  tillerSkeletonDir directory not exists, creating it ...')
                try:
                    os.mkdir(tillerSkeletonDir)
                except:
                    raise Exception('Cannot create tillerSkeletonDir directory: %s'%tillerSkeletonDir)
            fn_temp=block_i.replace('*','-')
            fnOutPrefix=os.path.join(tillerSkeletonDir,'TillerSkeleton_'+fn_temp)
            showOn=0
            TillerSkeletonPlot(plant_KP_virtualCoor_dict[block_i],fnOutPrefix,showOn) # plot tiller 2D skeleton for checking
        M_tiller_all=np.zeros([0,17])
        cultivar_i,plant_i_str,tiller_i_str=block_i.split('*')
        plant_i=int(plant_i_str)
        tiller_i=(plant_i-1)*10+int(tiller_i_str)
        stem_PTs=np.array(plant_KP_virtualCoor_dict[block_i]['stem'][0])*plantH_coef
        stem_PTs_num=len(stem_PTs)
        M_stem_all=np.zeros([0,17])
        barcode=[plant_i,tiller_i,0,3,0]
        stem_seg_len=[np.linalg.norm(np.array(stem_PTs[pt_i])-stem_PTs[pt_i-1]) for pt_i in range(1,stem_PTs_num-1)]
        stem_height_tot=sum(stem_seg_len)
        stem_node_diameter=[]
        for pt_i in range(0,stem_PTs_num-1):
            dist1=max(sum(stem_seg_len[0:pt_i]),0.00001)
            coef1=1/dist1
            coef2=1/max(abs(stem_height_tot/2-dist1),0.00001)
            coef3=1/max(stem_height_tot-dist1,0.00001)
            coefs=[coef1,coef2,coef3]
            node_diameter_temp=np.dot(stem_diameter,coefs)/sum([coef1,coef2,coef3])
            stem_node_diameter.append(node_diameter_temp)
        for pt_i in range(1,stem_PTs_num-1): # stem_PTs_num-1: do not use the last point (extendTip)
            ave_height=(sum(stem_seg_len[0:pt_i-1])+sum(stem_seg_len[0:pt_i]))/2
            if ave_height>stem_midH*plantH_coef:
                LNC_Kt_Kr=[stem_LNC_distribution[0],Kt_stem_distribution[0],Kr_stem_distribution[0]]
            else:
                LNC_Kt_Kr=[stem_LNC_distribution[1],Kt_stem_distribution[1],Kr_stem_distribution[1]]
            organ_base=[stem_PTs[pt_i-1][0],0,stem_PTs[pt_i-1][1]]
            organ_tip=[stem_PTs[pt_i][0],0,stem_PTs[pt_i][1]]
            stem_seg_diameter=[stem_node_diameter[pt_i-1],stem_node_diameter[pt_i]]
            M_stem = buildStem(barcode,LNC_Kt_Kr,organ_base,organ_tip,stem_seg_diameter) # plant, tiller, Position, organCode
            M_stem_all=np.vstack((M_stem_all,M_stem))
        M_tiller_all=np.vstack((M_tiller_all, M_stem_all))
        num_leaf=len(plant_KP_virtualCoor_dict[block_i]['leaf'])
        for leaf_i in range(num_leaf):
            leaf_order=num_leaf-leaf_i
            if leaf_order==1 and len(plant_KP_virtualCoor_dict[block_i]['ear'][0])>0: # ear presents, flag leaf
                widLen_relation_coef = widLen_relation_flag_coef
                widFunc = np.array(widFunc_flag)
            else:
                widLen_relation_coef = widLen_relation_other_coef
                widFunc = np.array(widFunc_other)
            barcode=[plant_i,tiller_i,leaf_order,4,0]
            LNC_Kt_Kr=[leaf_LNC_distribution[leaf_order-1],Kt_leaf_distribution[leaf_order-1],Kr_leaf_distribution[leaf_order-1]]
            leaf_PTs=np.array(plant_KP_virtualCoor_dict[block_i]['leaf'][leaf_i])
            leaf_base=np.array(leaf_PTs[0])*plantH_coef
            leaf_PTs=leaf_PTs-leaf_PTs[0]+leaf_base
            organ_PT_coords=[[leaf_PTs[i][0],0,leaf_PTs[i][1]] for i in range(len(leaf_PTs))]
            organ_PT_coords=np.array(organ_PT_coords)
            
            minDist_pos=0
            minDist=100000000
            for stemPt_i in range(stem_PTs_num):
                dist_temp=np.linalg.norm(leaf_base-stem_PTs[stemPt_i])
                if dist_temp<minDist:
                    minDist=dist_temp
                    minDist_pos=stemPt_i
            minDist_prev=max(0,minDist_pos-1)
            minDist_post=min(stem_PTs_num-1,minDist_pos+1)
            vec0=stem_PTs[minDist_pos] - leaf_base
            vec0_norm=np.linalg.norm(vec0)
            if abs(vec0_norm)<0.0001:
                stem_seg_base=[stem_PTs[minDist_pos][0],0,stem_PTs[minDist_pos][1]]
                stem_seg_tip=[stem_PTs[minDist_post][0],0,stem_PTs[minDist_post][1]]
            else:
                vec0=vec0/vec0_norm
                vec1=stem_PTs[minDist_prev] - leaf_base
                vec1=vec1/np.linalg.norm(vec1)
                vec2=stem_PTs[minDist_post] - leaf_base
                vec2=vec2/np.linalg.norm(vec2)
                val1=math.acos(max(-1,min(1,np.dot(vec0,vec1))))
                val2=math.acos(max(-1,min(1,np.dot(vec0,vec2))))
                if val1<val2: # at the same side of leaf line
                    stem_seg_base=[stem_PTs[minDist_pos][0],0,stem_PTs[minDist_pos][1]]
                    stem_seg_tip=[stem_PTs[minDist_post][0],0,stem_PTs[minDist_post][1]]
                else:
                    stem_seg_base=[stem_PTs[minDist_prev][0],0,stem_PTs[minDist_prev][1]]
                    stem_seg_tip=[stem_PTs[minDist_pos][0],0,stem_PTs[minDist_pos][1]]
            M_leaf = buildLeaf(barcode,LNC_Kt_Kr,organ_PT_coords,stem_seg_base,stem_seg_tip,widLen_relation_coef,widFunc,rollRad,sheathRollRad,transLen,ptDense,halfSplit,randZRot_max,tortAng,tortAng_start_pos,tortAng_end_pos,leafAng_coef)
            M_tiller_all=np.vstack((M_tiller_all,M_leaf))
        #print(np.array(plant_KP_virtualCoor_dict[block_i]['ear'][0]))
        #print(np.array(plant_KP_virtualCoor_dict[block_i]['ear'][0][0]))
        #raise Exception(np.array(plant_KP_virtualCoor_dict[block_i]['ear'][0][0])*plantH_coef)
        ear_PTs=np.array(plant_KP_virtualCoor_dict[block_i]['ear'][0])
        if len(ear_PTs) and withSpike:
            ear_PTs=ear_PTs[0]*plantH_coef + (ear_PTs-ear_PTs[0])*earLen_coef
            M_ear_all=np.zeros([0,17])
            barcode=[plant_i,tiller_i,-2,5,0]
            LNC_Kt_Kr=[rachis_LNC,Kt_rachis,Kr_rachis,spikelet_LNC,Kt_spikelet,Kr_spikelet,awn_LNC,Kt_awn,Kr_awn]
            organ_PT_coords=[[ear_PTs[i][0],0,ear_PTs[i][1]] for i in range(len(ear_PTs))]
            organ_PT_coords=np.array(organ_PT_coords)
            ear_PTs=np.array(ear_PTs)
            ear_len=sum([np.linalg.norm([ear_PTs[i]-ear_PTs[i+1]]) for i in range(len(ear_PTs)-1)])
            spikeletNum=round(float(ear_len*spikeletDense))
            M_ear = buildEar(barcode,LNC_Kt_Kr,organ_PT_coords,spikeletNum,spikeletBranchLen,spikeletBranchThick,awnLen,awnThick,grainLen,grainWid,grainThick,grainLen2,grainWid2,grainThick2,grainShrinkCoef,awnAng_coef_bot,awnAng_coef_mid,awnAng_coef_top,glumeAng,LEMMA_ANG,SPIKELET_ANG,zRotAng,earRachis_diameter,allowedLayerNum,glumeAwnLen)
            M_tiller_all=np.vstack((M_tiller_all, M_ear))
        M_multi_tillers.append(M_tiller_all)
    
    ### reconstruct the first plant
    M_plant=M_multi_tillers[0]
    print('  Writing Triangle-Patch into file...')
    format_str='%d\t'*5+'%.3f\t'*12
    format_str=format_str[0:-1]
    np.savetxt(fnOut_plant,M_plant,fmt=format_str) # , delimiter='\t'
    picOut_plant=fnOut_plant[0:-4]+'.png'
    command='python .\plantDraw3D.py -i %s -o %s -ccs 4 -cql 0 -cqu 1 -cm plasma -alpha 1 -va1 0 -va2 90'%(fnOut_plant,picOut_plant)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)
    ### reconstruct a canopy
    M_canopy = buildCanopy(M_multi_tillers,row_dist,sowing_width,plant_density,rowNum_simulate,rowLen_simulate,tillerAng,Kr_soil)
    print('  Writing Triangle-Patch into file...')
    format_str='%d\t'*5+'%.3f\t'*12
    format_str=format_str[0:-1]
    np.savetxt(fnOut,M_canopy,fmt=format_str) # , delimiter='\t'
    picOut_canopy=fnOut[0:-4]+'.png'
    command='python .\plantDraw3D.py -i %s -o %s -ccs 4 -cql 0 -cqu 1 -cm plasma -alpha 1 -va1 30 -va2 30'%(fnOut,picOut_canopy)
    try:
        subprocess.run(command, check = True)
    except:
        raise Exception('Command execution error: %s'%command)
    
    eclipse=time.time()-start
    print('  Canopy construction done! Time used: ', eclipse)