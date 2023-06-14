
#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Geophysical Research Letters 06/2023
# Author script: A. Boot
# Contact person: Amber Boot (she/they; d.boot@uu.nl)

# Script for selecting integration region
# Data necessary to run script can be downloaded from https://esgf-node.llnl.gov/search/cmip6/

#%% Import modules
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import xesmf as xe

#%% Directory where data is stored
datadir='/Users/daan/Desktop/Cimatoribus_model_amber/Repository/CMIP6 fit/Fit data/' 

#%% Data dir
datadir='/Users/daan/Desktop/ESM_data/' 
  
#%% Function to load in data and regrid onto 1x1 grid
def load_data(string1,string2,var1,var2):
    load_var1 = xr.open_dataset(f'{datadir}/'+string1)
    VAR1=load_var1[var1].compute().squeeze()
   
    grid1 = VAR1
    regridder1 = xe.Regridder(grid1, ds_out, 'bilinear',periodic=True,ignore_degenerate=True) 
    
    VR1 = regridder1(VAR1)
    
    load_var1 = xr.open_dataset(f'{datadir}/'+string2)
    VAR2=load_var1[var2].compute().squeeze()
    
    grid2 = VAR2
    regridder2 = xe.Regridder(grid2, ds_out, 'bilinear',periodic=True,ignore_degenerate=True) 
    
    VR2 = regridder2(VAR2)
    
    return VR1,VR2

#%% Function to load in data and regrid onto 1x1 grid
def load_data2(string1,var1):
    load_var1 = xr.open_dataset(f'{datadir}/'+string1)
    VAR1=load_var1[var1].compute().squeeze()
    
    grid1 = VAR1
    regridder1 = xe.Regridder(grid1, ds_out, 'bilinear',periodic=True,ignore_degenerate=True) 
    
    VR1 = regridder1(VAR1)
    
    return VR1

#%% Load data of a specific model
model = 'NorCPM1'
var1 = 'vsf'
var2 = 'evs'

mem = 1

string1 = var1+'_1pct_co2_'+model+'_'+str(mem)+'.nc'
string2 = var2+'_1pct_co2_'+model+'_'+str(mem)+'.nc'

ds_out = xe.util.grid_global(1, 1)
X = load_data2(string1,var1)

#%% Lat lon data
lat1=10
lat2=50
lat3=50+90
lat4=70+90
lat5=60
lat6=50+90

lon1=-81-180
lon2=-1-160

#%% Area data
data1='/Users/daan/CESM2_data'   # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
area_gr=load_var1['areacello'][:,:].compute().squeeze()

lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%% Mask representing box t in Box model
load_var1 = xr.open_dataset(f'{datadir}/thermocline_mask.nc')
mask=load_var1['vsf'][:,:].compute().squeeze()

#%% areas
area_s=np.array(area_gr[lat1:lat2,:])
area_n=np.array(area_gr[lat3:lat4,lon1:lon2])
area_t=np.array(area_gr[lat5:lat6,lon1:lon2])*np.array(mask)

#%% Function to integrate over box n, s and t
def integrate(data1,data2):
    VS1_s=data1[:,lat1:lat2,:]
    VS1_n=data1[:,lat3:lat4,lon1:lon2]
    VS1_t=data1[:,lat5:lat6,lon1:lon2]*np.array(mask)
    
    VS1_si=np.nanmean(np.nanmean(VS1_s*area_s,axis=1),axis=1)
    VS1_ni=np.nanmean(np.nanmean(VS1_n*area_n,axis=1),axis=1)
    VS1_ti=np.nanmean(np.nanmean(VS1_t*area_t,axis=1),axis=1)
    
    VS2_s=data2[:,lat1:lat2,:]
    VS2_n=data2[:,lat3:lat4,lon1:lon2]
    VS2_t=data2[:,lat5:lat6,lon1:lon2]*np.array(mask)

    VS2_si=np.nanmean(np.nanmean(VS2_s*area_s,axis=1),axis=1)
    VS2_ni=np.nanmean(np.nanmean(VS2_n*area_n,axis=1),axis=1)
    VS2_ti=np.nanmean(np.nanmean(VS2_t*area_t,axis=1),axis=1)
    
    VS3_ti = VS1_ti-VS2_ti
    VS3_si = VS1_si-VS2_si
    VS3_ni = VS1_ni-VS2_ni

    return VS3_ti,VS3_si,VS3_ni

#%% Function to integrate over box n, s and t
def integrate2(data1):
    VS1_s=data1[:,lat1:lat2,:]
    VS1_n=data1[:,lat3:lat4,lon1:lon2]
    VS1_t=data1[:,lat5:lat6,lon1:lon2]*np.array(mask)
    
    VS1_si=np.nansum(np.nansum(VS1_s*area_s,axis=1),axis=1)
    VS1_ni=np.nansum(np.nansum(VS1_n*area_n,axis=1),axis=1)
    VS1_ti=np.nansum(np.nansum(VS1_t*area_t,axis=1),axis=1)

    return VS1_ti,VS1_si,VS1_ni

#%% Integrate over region
tm = 60
P_1 = integrate2(X)

#%% Turn into xarray
data_xr = xr.DataArray(np.array(P_1), 
coords={'location': [1,2,3],'time': X['time']}, 
dims=["location", "time"])

#%% Save integrated file
filename = 'wfo_1pct_co2_'+model+'_'+str(mem)+'.nc'
data_xr.to_netcdf(datadir+filename)
