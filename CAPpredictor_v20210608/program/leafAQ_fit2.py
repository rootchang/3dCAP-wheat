#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
leafAQ_fit2.py

Fit leaf photosynthetic parameters (at 25 oC) based on measured AQ data by Licor-6400
"""

import argparse
import os
import xlrd
import time
from datetime import datetime
from xlrd import xldate_as_tuple
from datetime import timedelta
#import xlwt
from openpyxl import Workbook
import re
import numpy as np
import random
from scipy.interpolate import PchipInterpolator
from scipy.optimize import leastsq,curve_fit,fmin_slsqp
import math

global epsilon
epsilon=0.0001

__author__="CHANG Tiangen"
__date__="20201207"
__copyright__="Copyright 2020 Aspar CHANG"
__license__="SIPPE"
__version__="0.2"
__email__="changtiangen@sippe.ac.cn"
__status__="Prototype"
__credits__=[]

def createHelp():
    """
    Create the command line interface of the programme.
    """
  
    epilog_string="Any bug is welcome reported to changtiangen@sippe.ac.cn"
    description_string='The program is used to fit leaf AQ curve.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--di', dest='dir_in', help='Directory for input file(s) containing measured AQ data')
    parser.add_argument('-i2', '--paraFn', dest='para_fn', help='Input file containing Amax_net and Rd temperature response coefficients')
    parser.add_argument('-o', '--do', dest='dir_out', help='Directory for output file(s) containing fitted AQ data and parameters')
    op=parser.parse_args()
    return op

def uniqueFitJudger(itera_test_in,fit_error_new,fit_error,popt,popt_old,uniq_old):
    itera_test_out=itera_test_in
    fit_error_out=fit_error
    popt_new=popt_old
    #print(popt_new,popt_old)
    theta,Amax_net25,phi_CO2,Rd25=popt_old
    theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq=uniq_old
    if fit_error_new<fit_error-epsilon:
        #print('fit_error_new: %.2f'%fit_error_new)
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
        if abs(Amax_net25-popt[1])/(Amax_net25+epsilon)>0.05:
            Amax25_uniq='NO'
        if abs(phi_CO2-popt[2])/(phi_CO2+epsilon)>0.05:
            phi_CO2_uniq='NO'
        if abs(Rd25-popt[3])/(Rd25+epsilon)>0.05:
            Rd25_uniq='NO'
    uniq_new=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq]
    return itera_test_out,fit_error_out,popt_new,uniq_new

def fit_AQ(PAR_ave,A_ave,CT_ave,Amax_net_coef,Rd_coef):
    x_tot=np.array(PAR_ave)
    y=np.array(A_ave)
    z=np.array(CT_ave)
    ################ simple total absorbed PPFD to A AQ fitting ################
    itera_test=0
    fit_error=100000000
    theta=0
    Amax_net25=0
    phi_CO2=0
    Rd25=0
    theta_uniq='YES'
    Amax25_uniq='YES'
    phi_CO2_uniq='YES'
    Rd25_uniq='YES'
    
    while itera_test < 50:
        #print('            Simple AQ fitting: iteration %d, fit_error %.1f'%(itera_test,fit_error))
        theta_ini=random.uniform(0.01,0.999)
        Amax_net25_ini=random.randint(1,100)*1.01
        phi_CO2_ini=random.uniform(0,0.1)
        Rd25_ini=random.randint(1,100)*1.01
        paras0=np.array([theta_ini,Amax_net25_ini,phi_CO2_ini,Rd25_ini]) # theta,Amax_net,alpha,Rd
        popt= fmin_slsqp(fit_residuals, paras0, bounds=[(0.01,0.999),(0,101),(0,0.1),(0,101)], args=(x_tot,y,z,Amax_net_coef,Rd_coef),iprint = 0)
        popt=popt[0:]
        if math.isnan(popt[0]):
            continue
        if theta+Amax_net25+phi_CO2==0:
            theta,Amax_net25,phi_CO2,Rd25=popt
            fit_error=fit_residuals(popt[0:],x_tot,y,z,Amax_net_coef,Rd_coef)
        else:
            fit_error_new=fit_residuals(popt[0:],x_tot,y,z,Amax_net_coef,Rd_coef)
            popt_old=[theta,Amax_net25,phi_CO2,Rd25]
            uniq_old=[theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq]
            itera_test,fit_error,popt_new,uniq_new = uniqueFitJudger(itera_test,fit_error_new,fit_error,popt,popt_old,uniq_old)
            theta,Amax_net25,phi_CO2,Rd25=popt_new
            theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq=uniq_new
        itera_test=itera_test+1
    Rd_array=Rd25*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Amax_net_array=Amax_net25*(Amax_net_coef[0]*z**2 + Amax_net_coef[1]*z + Amax_net_coef[2]) # Rd-T empirical relation
    Ag_max=Amax_net_array+Rd_array
    A_fit_array=(phi_CO2*x_tot+Ag_max-np.sqrt((phi_CO2*x_tot+Ag_max)**2-4*theta*phi_CO2*x_tot*Ag_max))/(2*theta)-Rd_array
    
    R2_fit=1-sum((A_fit_array-y)**2)/sum((y-np.mean(y))**2)
    
    return theta,Amax_net25,phi_CO2,Rd25,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,fit_error,R2_fit,A_fit_array


def fit_residuals(paras,x_tot,y,z,Amax_net_coef,Rd_coef):#y=A_entire, x=PPFD, z=CT, y2=A_noAwn, z2=CT2
    theta,Amax_net25,phi_CO2,Rd25=paras
    Rd_facet_array=Rd25*(Rd_coef[0]*z**2 + Rd_coef[1]*z + Rd_coef[2]) # Rd-T empirical relation
    Amax_net_facet_array=Amax_net25*(Amax_net_coef[0]*z**2 + Amax_net_coef[1]*z + Amax_net_coef[2]) # A-T empirical relation
    Ag_facet_max=Amax_net_facet_array+Rd_facet_array
    
    sqrt1=(phi_CO2*x_tot+Ag_facet_max)**2-4*theta*phi_CO2*x_tot*Ag_facet_max
    y_fit=(phi_CO2*x_tot+Ag_facet_max-np.sqrt(sqrt1))/(2*theta)-Rd_facet_array
    fit_err=sum((y_fit-y)**2)
    #print('Simple fit: ',y_fit,y_entire)
    return fit_err
    

if __name__=="__main__":
    
    op=createHelp()
    
    photosynthesis_folder=op.dir_in
    photosynthesis_file_list=os.listdir(photosynthesis_folder)
    
    result_dir=op.dir_out
    config_Tresponse_fn=op.para_fn
    
    for line in open(config_Tresponse_fn,'r'):
        words=line.strip().split('\t')
        if words[0]=='coef_leaf_Asat_airT':
            Amax_net_coef=[float(c) for c in words[1:]]
        elif words[0]=='coef_leaf_Rd_airT':
            Rd_coef=[float(c) for c in words[1:]]
    
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    for fn in photosynthesis_file_list:
        #workbook = xlwt.Workbook(encoding='utf-8')
        workbook = Workbook()
        #worksheet = workbook.add_sheet('1')
        worksheet = workbook.active
        record=['Remark','PARi','airT','blockT','leafT','CO2R','A','A_fit','Amax_gross25','theta','alpha','Rd25','Amax_net25','AmaxUniq','thetaUniq','alphaUniq','RdUniq','R2','err']
        worksheet.append(record)
        write_line=1
        
        print(fn+' in processing...')
        data=xlrd.open_workbook(photosynthesis_folder+'\\'+fn)
        Num_sheet=len(data.sheet_names())
        for sheet_i in range(Num_sheet):
            table=data.sheet_by_index(sheet_i)
            row_num=table.nrows
            PARi_list=[]
            CO2R_list=[]
            AT_list=[]
            BT_list=[]
            LT_list=[]
            A_list=[]
            NA_old=''
            for i in range(row_num):
                words=table.row_values(i)
                if words[0]:
                    NA=words[0]
                    if len(A_list):
                        theta,Amax_net25,phi_CO2,Rd25,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,fit_error,R2_fit,Afit=fit_AQ(PARi_list,A_list,BT_list,Amax_net_coef,Rd_coef)
                        Amax25_gross=Amax_net25+Rd25
                        worksheet.append([NA_old,PARi_list[0],AT_list[0],BT_list[0],LT_list[0],CO2R_list[0],A_list[0],Afit[0],Amax25_gross,theta,phi_CO2,Rd25,Amax_net25,Amax25_uniq,theta_uniq,phi_CO2_uniq,Rd25_uniq,R2_fit,fit_error])
                        for item_i in range(1,len(A_list)):
                            worksheet.append(['',PARi_list[item_i],AT_list[item_i],BT_list[item_i],LT_list[item_i],CO2R_list[item_i],A_list[item_i],Afit[item_i]])
                        PARi_list=[]
                        CO2R_list=[]
                        AT_list=[]
                        BT_list=[]
                        LT_list=[]
                        A_list=[]
                    NA_old=NA
                try:
                    float(words[1])
                except:
                    continue
                #raise Exception([float(c) for c in [words[19]]+words[9:13]+[words[1]]])
                PARi,AT,BT,LT,CO2R,A=[float(c) for c in [words[19]]+words[9:13]+[words[1]]]
                if len(PARi_list) and PARi>PARi_list[-1]+20:
                    PARi_list[-1]=PARi
                    AT_list[-1]=AT
                    BT_list[-1]=BT
                    LT_list[-1]=LT
                    CO2R_list[-1]=CO2R
                    A_list[-1]=A
                else:
                    PARi_list.append(PARi)
                    AT_list.append(AT)
                    BT_list.append(BT)
                    LT_list.append(LT)
                    CO2R_list.append(CO2R)
                    A_list.append(A)
                if i==row_num-1: # last line
                    if len(A_list):
                        theta,Amax_net25,phi_CO2,Rd25,theta_uniq,Amax25_uniq,phi_CO2_uniq,Rd25_uniq,fit_error,R2_fit,Afit=fit_AQ(PARi_list,A_list,BT_list,Amax_net_coef,Rd_coef)
                        Amax25_gross=Amax_net25+Rd25
                        worksheet.append([NA_old,PARi_list[0],AT_list[0],BT_list[0],LT_list[0],CO2R_list[0],A_list[0],Afit[0],Amax25_gross,theta,phi_CO2,Rd25,Amax_net25,Amax25_uniq,theta_uniq,phi_CO2_uniq,Rd25_uniq,R2_fit,fit_error])
                        for item_i in range(1,len(A_list)):
                            worksheet.append(['',PARi_list[item_i],AT_list[item_i],BT_list[item_i],LT_list[item_i],CO2R_list[item_i],A_list[item_i],Afit[item_i]])

        workbook.save(os.path.join(result_dir,fn[0:-5]+'_out.xlsx'))
