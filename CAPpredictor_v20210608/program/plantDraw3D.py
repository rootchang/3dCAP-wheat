#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plantDraw3D.py

draw 3D plant/canopy from triangulation data file
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.colors as colors
import argparse
from pylab import get_cmap
from matplotlib.tri import Triangulation, LinearTriInterpolator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib
import time
import json
#import cStringIO
from matplotlib.colors import LightSource

__author__="CHANG Tiangen"
__date__="20200918"
__copyright__="Copyright 2020 Aspar CHANG"
__license__="CEMPS-CAS"
__version__="0.1"
__email__="changtiangen@cemps.ac.cn"
__status__="Prototype"
__credits__=[]


def createHelp():
    """
    Create the command line interface of the programme.
    """
    
    epilog_string="Any bug is welcome reported to changtiangen@cemps.ac.cn"
    description_string='The program is going to draw a given 3D plant/canopy.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--input file', dest='fnIn', type=str, help='<Required> input triangulation data file', required=True)
    parser.add_argument('-o', '--output file', dest='fnOut', type=str, help='<Required> output 3D image (.png)', required=True)
    parser.add_argument('-ccs', '--color column', dest='ccs', nargs='*', type=int, default=[4],help='Coloring facet based on values of [col_i1,col_i2,...], respectively, e.g., [4] (based on fnIn column4 "organCode"), [-1] (based on fnIn last column magnitude), [0] (default jet)')
    parser.add_argument('-cm', '--color mode', dest='cm', type=str, default='plasma',help='Color map code, e.g., plasma,inferno,winter,summer,jet')
    parser.add_argument('-axLim', '--axes range', dest='axLim', nargs='*',type=float, default=[], help='x,y,z axes range')
    parser.add_argument('-alpha', '--transparence', dest='alpha', type=float,default=1, help='Transparence of patches')
    parser.add_argument('-va1', '--view angles1', dest='va1', nargs='*',type=int,default=[45], help='View angle 1. Rotating angles around y axis')
    parser.add_argument('-va2', '--view angles2', dest='va2', nargs='*',type=int,default=[45], help='View angle 2. Rotating angles around z axis')
    parser.add_argument('-cql', '--color quantileLow', dest='cql',type=float,default=0.05, help='Lower quantile for coloring')
    parser.add_argument('-cqu', '--view quantileUpper', dest='cqu',type=float,default=0.95, help='Upper quantile for coloring')
    parser.add_argument('-header', '--header', dest='header',type=int,default=0, help='The input file has a header line (1) or not (0)')
    parser.add_argument('-show', '--show', dest='showOn',type=int,default=0, help='To show the draw on screen (1) or save it to file (0)')

    op=parser.parse_args()
    return op
    
def fileLineNumFastCounter(file_name):
    from itertools import (takewhile, repeat)
    buffer = 1024 * 1024
    with open(file_name) as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        return sum(buf.count('\n') for buf in buf_gen)

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
        
def plot3Dpatch(M_plant,colorCols,colorMode,ax_range,t_alpha,viewAngles1,viewAngles2,color_ql,color_qu,fnOut,showOn):
    # fonts setting
    ticklabelFont=10
    labelFont=12
    legendFont=10
    refFont=10
    
    triangle_vertices=np.array([[c[5:8],c[8:11],c[11:14]] for c in M_plant])
    col_num=len(M_plant[0,:])
    if not ax_range:
        x_all=np.hstack((M_plant[:,5],M_plant[:,8],M_plant[:,11]))
        y_all=np.hstack((M_plant[:,6],M_plant[:,9],M_plant[:,12]))
        z_all=np.hstack((M_plant[:,7],M_plant[:,10],M_plant[:,13]))
        xmin=min(x_all)
        ymin=min(y_all)
        zmin=max(0,min(z_all))
        xmax=max(x_all)
        ymax=max(y_all)
        zmax=max(z_all)*1.1
    else:
        [xmin,xmax,ymin,ymax,zmin,zmax]=ax_range
    
    for colorCol in colorCols:
        if colorCol<1:
            colorCol=col_num+colorCol+1
        fig = plt.figure(figsize=(6,6))
        left_edge=0.
        bottom_edge=0.
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
            ax_list.append(plt.axes(rect_line[i], projection='3d'))
            ax_list[i].tick_params(direction='in',axis='both',length=2, which='major', labelsize=ticklabelFont,pad=2,top='off', right='off')
            ax_list[i].yaxis.set_label_coords(-0.1, 0.5)
            ax_list[i].xaxis.set_label_coords(0.5, -0.1)
    
        ax_list[0].grid(False)
        ax_list[0].set_axis_off()
        color_list=M_plant[:,colorCol-1]
        
        # black the bg
#         bg=M_plant[:,2]>-6
#         color_list=color_list*bg
        
        color_totNum=len(set(color_list))
        #print(color_totNum)
        if color_totNum==1: # plant_id
            color_list[-1]=color_qu
            color_totNum=1000
        cmin=0
        cmax=1
        norm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
        if colorCol==4: # organType id
            colormap=[[128/255,128/255,128/255],[50/255,205/255,50/255], [50/255,205/255,50/255],[0/255,100/255,0/255], [255/255,185/255,15/255]] # 50/255,205/255,50/255
            patchColor_mat=[colormap[int(c)-1]+[1] for c in color_list]
            cmap_name='organ_cmap'
            cm = colors.LinearSegmentedColormap.from_list(cmap_name, colormap, N=5)
            bounds = [ round(elem, 2) for elem in np.linspace(cmin, cmax, 5)] # 
            tick_list=[ round(elem, 2) for elem in np.linspace(cmin+0.125, cmax-0.125, 4)] # 
            label_list=['Soil','culm','Leaf','Spike']
        else:
            cm = get_cmap(colorMode)
            if color_qu<=1:
                color_pmin=np.quantile(color_list,color_ql)
                color_pmax=np.quantile(color_list,color_qu)
            else:
                color_pmin=np.array(color_ql)
                color_pmax=np.array(color_qu)
            color_list_trunc=[min(max(c,color_pmin),color_pmax) for c in color_list]
            if color_pmax>color_pmin:
                color_list_norm=(color_list_trunc-color_pmin)/(color_pmax-color_pmin) # norm to [0,1]
            else:
                color_list_norm=(color_list_trunc-color_pmin)
            #print(color_pmin,color_pmax,min(color_list),max(color_list))
            patchColor_mat=cm(color_list_norm)
            #print(set(color_list_norm))
            bounds = [ round(elem, 2) for elem in np.linspace(cmin, cmax, color_totNum)] # 
            tick_list=[ round(elem, 2) for elem in np.linspace(cmin, cmax, 6)] # 
            label_list=np.array(tick_list)*(color_pmax-color_pmin)+color_pmin
            if label_list[-1]<10:
                label_list=[str(round(c,1)) for c in label_list]
            else:
                label_list=[str(int(round(c,0))) for c in label_list]
            label_list[0]='<'+label_list[0]
            label_list[-1]='>'+label_list[-1]
        collection = Poly3DCollection(triangle_vertices,facecolors=patchColor_mat,edgecolors=patchColor_mat, linewidth=0.2,antialiased = False, alpha=t_alpha) #  norm=norm,
        #collection = Poly3DCollection(triangle_vertices,facecolors='k',edgecolors='None', linewidth=0,antialiased = True, alpha=0.6) #  norm=norm,
#         collection = Poly3DCollection(triangle_vertices,facecolors='k',edgecolors='k', linewidth=0.2,antialiased = True, alpha=1) #  norm=norm,
        #ls = LightSource(azdeg=225.0, altdeg=45.0)
        p=ax_list[0].add_collection(collection)
        ax3 = fig.add_axes([0.85, 0.3, 0.02, 0.5]) # left, bottom, width, height
        #  to use 'extend', you must specify two extra boundaries:
        #  boundaries= [1.2] + bounds + [2.6]
        boundaries=bounds
        cb3 = matplotlib.colorbar.ColorbarBase(ax3, cmap=cm,
                                     boundaries=bounds,
                                     ticks=tick_list, # optional
                                     spacing='proportional',
                                     orientation='vertical')
        ax3.set_yticklabels(label_list)
        ax_list[0].set_xlim([xmin,xmax])
        ax_list[0].set_ylim([ymin,ymax])
        ax_list[0].set_zlim([zmin,zmax])
        axisEqual3D(ax_list[0])
        # rotate the axes and update
        for angle1 in viewAngles1:
            for angle2 in viewAngles2:
                print('Color column %d, angles %d-%d in processsing ...'%(colorCol,angle1,angle2))
                ax_list[0].view_init(angle1, angle2)
                #plt.draw()
                #plt.pause(.001)
                if showOn:
                    plt.show()
                else:
                    fnOut_temp=fnOut[0:-4]+'_'+'%02d'%colorCol+'-'+'%03d'%int(angle1)+'_'+'%03d'%int(angle2)+fnOut[-4:]
                    plt.savefig(fnOut_temp,format="png",dpi=600)
        plt.close()

if __name__=="__main__":
    start=time.time()
    op=createHelp()
    
    width_format=94
    welcome_info='Welcome to use plantDraw3D'
    time_info=time.asctime( time.localtime(time.time()) )
    print(welcome_info.center(width_format, '-'))
    print(time_info.center(width_format, ' '))
    
    fnIn=op.fnIn
    fnOut=op.fnOut
    colorCols=op.ccs
    colorMode=op.cm
    ax_range=op.axLim
    t_alpha=op.alpha
    viewAngle1=op.va1
    viewAngle2=op.va2
    color_ql=op.cql
    color_qu=op.cqu
    header=op.header
    showOn=op.showOn
    data = json.dumps(vars(op), indent=4,ensure_ascii=False, sort_keys=False,separators=(',', ':'))
    print('Running parameters:')
    print(data)
    
    if ax_range:
        if len(ax_range)!=6:
            raise Exception('axLim element number error: %d. Should be 6.'%len(ax_range))
    print('Processing input plant triangulation file...')
    triNum=fileLineNumFastCounter(fnIn)
    if header:
        triNum=triNum-1
    dataIn=open(fnIn,'r').readline()
    words=dataIn.strip().split('\t')
    colNum=len(words)
    dataIn=open(fnIn,'r')
    
    M_plant=np.zeros([triNum,colNum])
    i=0
    for line in dataIn:
        words=line.strip().split('\t')

        try:
            words=[float(c) for c in words]
            if 1:#not words[2]>-6:#words[24]>10 and words[24]<800:#1:#words[2]>-6: # and words[-1]>100:
                M_plant[i,:]=words
                i+=1
        except:
            continue
    M_plant=M_plant[0:i,:]
    print('Plotting the 3D patches ...')
    plot3Dpatch(M_plant,colorCols,colorMode,ax_range,t_alpha,viewAngle1,viewAngle2,color_ql,color_qu,fnOut,showOn)

    eclipse=time.time()-start
    print('All done! Time used: %.2fs.\n'%eclipse)