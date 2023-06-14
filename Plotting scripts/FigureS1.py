
#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Geophysical Research Letters 06/2023
# Author script: A. Boot
# Contact person: Amber Boot (she/they; d.boot@uu.nl)

# Script for plotting Figure S1

#%% Import modules
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#%% Data dir
datadir='/Users/daan/Desktop/Cimatoribus_model_amber/Repository/CMIP6 fit/Fit data/' 

#%% Whether to save figure and for what quality
save_fig = 'no' # 'yes' if figure needs to be saved
quality = 300 # Figure quality

#%% Function to load in data
def load_data(var,model,mem):
    load_var1 = xr.open_dataset(f'{datadir}/'+var+'_1pct_co2_'+model+'_'+str(mem)+'.nc')
    VAR1=load_var1[var1].compute().squeeze()
    
    return VAR1

#%% Load data based on wfo (water flux ocean)
var = 'wfo'
var1 = '__xarray_dataarray_variable__'

ACCESS_CM2_1 = load_data(var,'ACCESS-CM2',1)
ACCESS_ESM1_1 = load_data(var,'ACCESS-ESM1-5',1)

CanESM_1 = load_data(var,'CanESM5',32)

CanESM_1_11 = load_data(var,'CanESM5-1',11)
CanESM_1_12 = load_data(var,'CanESM5-1',12)
CanESM1_EM = (xr.concat([CanESM_1_11,CanESM_1_12],dim='member')).mean('member')

CMCC_HR_1 = load_data(var,'CMCC-CM2-HR4',1)
CMCC_SR_1 = load_data(var,'CMCC-CM2-SR5',1)
CMCC_ESM_1 = load_data(var,'CMCC-ESM2',1)

CNRM_CM_1 = load_data(var,'CNRM-CM6-1',12)
CNRM_HR_1 = load_data(var,'CNRM-CM6-1-HR',12)

CNRM_ESM_12 = load_data(var,'CNRM-ESM2-1',12)
CNRM_ESM_22 = load_data(var,'CNRM-ESM2-1',22)
CNRM_ESM_32 = load_data(var,'CNRM-ESM2-1',32)
CNRM_ESM_42 = load_data(var,'CNRM-ESM2-1',42)
CNRM_ESM_52 = load_data(var,'CNRM-ESM2-1',52)
CNRM_ESM_62 = load_data(var,'CNRM-ESM2-1',62)
CNRM_ESM_72 = load_data(var,'CNRM-ESM2-1',72)
CNRM_ESM_82 = load_data(var,'CNRM-ESM2-1',82)
CNRM_ESM_92 = load_data(var,'CNRM-ESM2-1',92)
CNRM_ESM_102 = load_data(var,'CNRM-ESM2-1',102)
CNRM_ESM_EM = (xr.concat([CNRM_ESM_12,CNRM_ESM_22,CNRM_ESM_32,CNRM_ESM_42,CNRM_ESM_52,CNRM_ESM_62,CNRM_ESM_72,CNRM_ESM_82,CNRM_ESM_92,CNRM_ESM_102],dim='member')).mean('member')

GFDL_CM_1 = load_data(var,'GFDL-CM4',1)
GFDL_ESM_1 = load_data(var,'GFDL-ESM4',1)

GISS_1G_1 = load_data(var,'GISS-E2-1-G',1)
GISS_1G_3 = load_data(var,'GISS-E2-1-G',3)
GISS_1G_5 = load_data(var,'GISS-E2-1-G',5)
GISS_1G_101 = load_data(var,'GISS-E2-1-G',101)
GISS_1G_102 = load_data(var,'GISS-E2-1-G',102)
GISS_1G_EM = (xr.concat([GISS_1G_1,GISS_1G_3,GISS_1G_5,GISS_1G_101,GISS_1G_102],dim='member')).mean('member')

GISS_2G_1 = load_data(var,'GISS-E2-2-G',1)

IPSL_5_1 = load_data(var,'IPSL-CM5A2-INCA',1)
IPSL_6_1 = load_data(var,'IPSL-CM6A-LR',1)

MIROC_ES2L_1 = load_data(var,'MIROC-ES2L',2)
MIROC_6_1 = load_data(var,'MIROC6',1)

MPI_HMA_1 = load_data(var,'MPI-ESM-1-2-HAM',1)

MRI_1 = load_data(var,'MRI-ESM2-0',1)
MRI_2 = load_data(var,'MRI-ESM2-0',2)
MRI_EM = (xr.concat([MRI_1,MRI_2],dim='member')).mean('member')

NESM_1 = load_data(var,'NESM3',1)

NorESM_LM_1 = load_data(var,'NorESM2-LM',1)
NorESM_MM_1 = load_data(var,'NorESM2-MM',1)

#%% Load data based on virtual salt fluxes
var = 'vsf'
var1 = '__xarray_dataarray_variable__'

CESM_1 = load_data(var,'CESM2',1)
CESM_W_1 = load_data(var,'CESM2-WACCM',1)
CESM_WF_1 = load_data(var,'CESM2-WACCM-FV2',1)

FGOALS_f_1 = load_data(var,'FGOALS-f3-L',1)
FGOALS_f_2 = load_data(var,'FGOALS-f3-L',2)
FGOALS_f_3 = load_data(var,'FGOALS-f3-L',3)
FGOALS_f_EM = (xr.concat([FGOALS_f_1,FGOALS_f_2,FGOALS_f_3],dim='member')).mean('member')

FGOALS_g_2 = load_data(var,'FGOALS-g3',2)
FGOALS_g_3 = load_data(var,'FGOALS-g3',3)
FGOALS_g_EM = (xr.concat([FGOALS_g_2,FGOALS_g_3],dim='member')).mean('member')

NorCPM_1 = load_data(var,'NorCPM1',1)     

#%% Store all model data in single area and get sign correct if necessary
ALL_MODELS = np.zeros((28,3,1800))     

ALL_MODELS[0,:,:] = ACCESS_CM2_1[:,:1800]
ALL_MODELS[1,:,:] = ACCESS_ESM1_1[:,:1800]
ALL_MODELS[2,:,:] = CanESM_1[:,:1800]
ALL_MODELS[3,:,:] = CanESM1_EM[:,:1800]
ALL_MODELS[4,:,:] = CESM_1[:,:1800]
ALL_MODELS[5,:,:] = CESM_W_1[:,:1800]
ALL_MODELS[6,:,:] = CESM_WF_1[:,:1800]
ALL_MODELS[7,:,:] = CMCC_HR_1[:,:1800]*-1
ALL_MODELS[8,:,:] = CMCC_SR_1[:,:1800]*-1

ALL_MODELS[9,:,:] = CMCC_ESM_1[:,:1800]*-1
ALL_MODELS[10,:,:] = CNRM_CM_1[:,:1800]*-1
ALL_MODELS[11,:,:] = CNRM_ESM_EM[:,:1800]*-1
ALL_MODELS[12,:,:] = FGOALS_f_EM[:,:1800]*-1
ALL_MODELS[13,:,:] = FGOALS_g_EM[:,:1800]*-1
ALL_MODELS[14,:,:] = GFDL_CM_1[:,:1800]
ALL_MODELS[15,:,:] = GFDL_ESM_1[:,:1800]
ALL_MODELS[16,:,:] = GISS_1G_EM[:,:1800]*-1
ALL_MODELS[17,:,:] = GISS_2G_1[:,:1800]*-1

ALL_MODELS[18,:,:] = IPSL_5_1[:,:1800]*-1
ALL_MODELS[19,:,:] = IPSL_6_1[:,:1800]*-1
ALL_MODELS[20,:,:] = MIROC_ES2L_1[:,:1800]
ALL_MODELS[21,:,:] = MIROC_6_1[:,:1800]*-1
ALL_MODELS[22,:,:] = MPI_HMA_1[:,:1800]
ALL_MODELS[23,:,:] = MRI_EM[:,:1800]
ALL_MODELS[24,:,:] = NESM_1[:,:1800]
ALL_MODELS[25,:,:] = NorCPM_1[:,:1800]*-35
ALL_MODELS[26,:,:] = NorESM_LM_1[:,:1800]

ALL_MODELS[27,:,:] = NorESM_MM_1[:,:1800]

#%% Determine multi model mean
MM_1 = np.mean(ALL_MODELS,axis=0)

STD_MM_1 = np.std(ALL_MODELS,axis=0)

MM_1_d1 = MM_1+STD_MM_1
MM_1_d2 = MM_1+2*STD_MM_1
MM_1_d3 = MM_1-STD_MM_1
MM_1_d4 = MM_1-2*STD_MM_1

#%% Define moving mean functions
tm=60 

# All models
def mov_mean1(data):
    dvs1 = pd.DataFrame(np.squeeze(data[:,0,:1800]))
    rolling_windows = dvs1.rolling(tm,axis=1)
    south=np.squeeze(rolling_windows.mean())
    
    dvs1 = pd.DataFrame(np.squeeze(data[:,1,:1800]))
    rolling_windows = dvs1.rolling(tm,axis=1)
    north=np.squeeze(rolling_windows.mean())
    
    return north, south   

# Multi model mean
def mov_mean(data):
    dvs1 = pd.DataFrame(data[0,:1800])
    rolling_windows = dvs1.rolling(tm)
    south=np.squeeze(rolling_windows.mean())
    
    dvs1 = pd.DataFrame(data[1,:1800])
    rolling_windows = dvs1.rolling(tm)
    north=np.squeeze(rolling_windows.mean())
    
    return north, south   

#%% Take moving mean
[am_north,am_south] = mov_mean1(-ALL_MODELS/1e9/2)
[mm_north,mm_south] = mov_mean(-MM_1/1e9/2)
[mm_d1_north,mm_d1_south] = mov_mean(-MM_1_d1/1e9/2)
[mm_d2_north,mm_d2_south] = mov_mean(-MM_1_d2/1e9/2)
[mm_d3_north,mm_d3_south] = mov_mean(-MM_1_d3/1e9/2)
[mm_d4_north,mm_d4_south] = mov_mean(-MM_1_d4/1e9/2)  

#%% Load in CO2 data
load_var1 = xr.open_dataset(f'{datadir}/co2mass_1pct_co2_CESM2_1.nc')
VAR1=load_var1['co2mass'].compute().squeeze()
CO2_1pct = VAR1/5.14e18*1.0e6 * 28.966 / 44.0

#%% Define fitting function
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

def func_lin(x,a,b):
    return a*x+b

def func_log(x,a,b):
    return a+b*np.log(x)

def fit(data):
    popt_1, pcov_1 = curve_fit(func_lin, np.array(CO2_1pct)[60:], np.array(data)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))
    popt_2, pcov_2 = curve_fit(func_log, np.array(CO2_1pct)[60:], np.array(data)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))

    fit1 = popt_1[0]*CO2_1pct+popt_1[1]
    fit2 = popt_2[0]+popt_2[1]*np.log(CO2_1pct)
    
    return fit1, fit2

#%% Define plot function
MS = 10
LW = 6
FS = 20

models = (['ACCESS-CM2-1','ACCESS-ESM1-5','CanESM5','CanESM5-1','CESM2','CESM2-WACCM','CESM2-WACCM-FV2','CMCC-CM2-HR4','CMCC-CM2-SR5','CMCC-ESM2','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','FGOALS-f3-L','FGOALS-g3','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','IPSL-CM5A2-INCA','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR','MRI-ESM2-0','NESM3','NorCPM1','NorESM2-LM','NorESM2-MM','Multimodel mean'])

def plot_line1(data,loc,i):
    plt.scatter(CO2_1pct,-data,s=MS,zorder=4,label='Data')
    [fit1,fit2]=fit(data)
    plt.plot(CO2_1pct,-fit1,color='tab:orange',linewidth=LW,label='Linear fit')
    plt.plot(CO2_1pct,-fit2,color='tab:green',linewidth=LW,label='Logaritmic fit')
    
    plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
    plt.ylabel('E$_s$ [Sv]',fontsize=FS-2)
    plt.title(models[i]+' '+loc,fontsize=FS)
    plt.xticks(fontsize=FS-3)
    plt.yticks(fontsize=FS-3)
    plt.legend(fontsize=FS-6)
    plt.grid()                                

#%% Create linear and log fit
data1=mm_north
popt_1, pcov_1 = curve_fit(func_lin, np.array(CO2_1pct)[60:], np.array(data1)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))
popt_2, pcov_2 = curve_fit(func_log, np.array(CO2_1pct)[60:], np.array(data1)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))

fit1 = popt_1[0]*CO2_1pct+popt_1[1]
fit2 = popt_2[0]+popt_2[1]*np.log(CO2_1pct)

r21=r2_score(np.array(mm_north)[60:], fit1[60:])
r22=r2_score(np.array(mm_north)[60:], fit2[60:])

print('Linear fit:')
print(str(np.round(popt_1[1],3))+' + ' + str(np.round(popt_1[0],6)) +' x CO2')

print('Log fit:')
print(str(np.round(popt_2[0],3))+' + ' + str(np.round(popt_2[1],3)) +' x ln(CO2)')

#%% Figure S1b
fig = plt.figure(figsize=(7, 5)) 
plt.scatter(CO2_1pct,-data1,s=MS,zorder=4,label='MM')
[fit1,fit2]=fit(data1)
plt.plot(CO2_1pct,-fit1,color='tab:orange',linewidth=LW,label='Linear fit; R$^2$ = ' +str(np.round(r21,2)))
plt.plot(CO2_1pct,-fit2,color='tab:green',linewidth=LW,label='Logarithmic fit; R$^2$ = ' +str(np.round(r22,2)))

plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
plt.ylabel('E$_s$ [Sv]',fontsize=FS-2)
plt.title('Multimodel mean',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-6)
plt.grid()

if save_fig == 'yes':
    plt.savefig('Figure_S1b.png', format='png', dpi=quality, bbox_inches = 'tight')
 
#%% Figure S1a
fig = plt.figure(figsize=(7, 5))
for i in range(28):
    plt.scatter(CO2_1pct,-np.array(am_north)[i,:],s=5)

plt.scatter(CO2_1pct,-np.array(mm_north),s=10,color='black',label=models[-1])
plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
plt.ylabel('E$_s$ [Sv]',fontsize=FS-2)
plt.title('CMIP6 ensemble',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-6,ncol=1,loc=2)
plt.grid()

if save_fig == 'yes':
    plt.savefig('Figure_S1a.png', format='png', dpi=quality, bbox_inches = 'tight')
    
#%% Figure S1c
CO2_new = np.arange(0,1001,10)
fit3 = (np.log((CO2_new+250)/300)*0.75+0.1)
fit2_new = 0.099 +0.079*np.log(CO2_new)
fit1_new = 0.00012*CO2_new+0.526

fig = plt.figure(figsize=(7, 5)) 

plt.plot(CO2_new,fit1_new,color='tab:orange',linewidth=LW,label='Eq. 1')
plt.plot(CO2_new,fit2_new,color='tab:green',linewidth=LW,label='Eq. 2')
plt.plot(CO2_new,fit3,color='tab:purple',linewidth=LW,label='Eq. 3')

plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
plt.ylabel('E$_s$ [Sv]',fontsize=FS-2)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-6)
plt.grid()
plt.xlim([0,1000])

if save_fig == 'yes':
    plt.savefig('Figure_S1c.png', format='png', dpi=quality, bbox_inches = 'tight')
 