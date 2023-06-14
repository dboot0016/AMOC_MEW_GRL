#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Geophysical Research Letters 06/2023
# Author script: A. Boot
# Contact person: Amber Boot (she/they; d.boot@uu.nl)

# Script for plotting Figure S3

#%% Import modules
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#%% Directory with data
Datadir='/Users/daan/Desktop/Cimatoribus_model_Amber/Repository/Data/'

#%% Plotting variables
varx = 'E$_a$' # Variable on x-axis plots
vary = 'CO$_2$' # Variable on y-axis plots; For AMOC plots: 'D', for CO2 plots: 'CO$_2$'
unitx = '[Sv]' # Unit on x-axis
vers = '2' # What Es-coupling equation to use: '2' or '3'

save_fig = 'no' # Whether to save figure: 'yes' to save
quality = 300 # Quality figure

LW = 4 # Linewidth of line
MS = 150 # Markersize
FS = 20 # Base fontsize in plots

#%% Model constants
D_dim = 1000 # Dimensionalize pycnocline depth
S_dim = 35 # Dimensionalize salinity

eta = 3e4 # Hydraulic constant
alpha = 2e-4 # Thermal expansion coefficient
beta = 8e-4 # Haline contraction coefficient
Tts = 10 # Temperature box ts
Tn = 5 # Temperature box n

tau = 0.1 # Wind stress
LxS = 3e7 # Zonal length Southern Ocean
rho0 = 1027.5 # Base density
fs = -1e-04 # Coriolis parameter
Agm = 1700 # Eddy diffusivity 
LxA = 1e7 # Zonal length Atlantic basin
Ly = 1e6 # Meridional extent frontal zone Southern Ocean

V0 = 3e17 # Total volume domain      
Vn = 3e15 # Volume box n       
Vs = 9e15 # Volume box s
A = 1e14 # Surface area box t 

epsilon = 1e-12 # Parameter for heavisede function
heaviside = lambda x: 0.5*(np.tanh(x/epsilon)+1) # Heaviside function

#%% Variables to load in data
a1 = np.repeat(np.arange(0,9.5,1),4)
a2 = np.tile([0,2,5,7],10)
a3 = np.tile([0,5],20)

#%% Create empty arrays for saddle node locations
L1 = np.zeros((41,))
L2 = np.zeros((41,))

#%% Select saddle node locations
for i in range(40):
    run1 = pd.read_csv(f'{Datadir}/'+'b.uncoupled_es_'+str(int(a1[i]))+str(a2[i])+str(a3[i]),sep='\s+',header=None)   
    run2 = pd.read_csv(f'{Datadir}/'+'b.uncoupled_es2_'+str(int(a1[i]))+str(a2[i])+str(a3[i]),sep='\s+',header=None)
    
    EA1 = run1[4]
    EA2 = run2[4]
    
    L1[i] = max(EA1)
    L2[i] = min(EA2)

#%% Saddle node location for last point
run1 = pd.read_csv(f'{Datadir}/'+'b.uncoupled_es_1000',sep='\s+',header=None)   
run2 = pd.read_csv(f'{Datadir}/'+'b.uncoupled_es2_1000',sep='\s+',header=None)

EA1 = run1[4]
EA2 = run2[4]

L1[-1] = max(EA1)
L2[-1] = min(EA2)

#%% Es-CO2 couplings
Es = np.arange(0,1.01,0.025)

CO2_1 = (Es-0.562)/(0.00012)
CO2_2 = np.exp(((Es-0.099)/0.069))
CO2_3 = np.exp((Es-0.1)/0.75)*300-250
 
#%% Figure S3
plt.figure(figsize=(8, 8))
fig, ax1 = plt.subplots(figsize=(8, 7.25))

plt.plot(Es,L1,'-.',color='tab:blue',linewidth=LW,label='On branch')
plt.plot(Es,L2,'--',color='tab:blue',linewidth=LW,label='Off branch')

plt.xlim([0,1.0])
plt.xlabel('E$_s$ [Sv]',fontsize=FS-3)
plt.ylabel('E$_a$ [Sv]',color='tab:blue',fontsize=FS-3)

plt.xticks([0,0.15,0.3,0.45,0.6,0.75,0.9],fontsize=FS-4)
plt.yticks([-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6],color='tab:blue',fontsize=FS-4)    
plt.grid()
plt.legend(fontsize=FS-6,loc=2)

ax2 = ax1.twinx()
ax2.plot(Es,CO2_3,linewidth=LW,color='g')
ax2.set_ylabel('CO$_2$ [ppm]', color='g',fontsize=FS-3)
plt.yticks([0,100,200,300,400,500,600,700],color='g',fontsize=FS-4)    

plt.xlim([0,1.0])

if save_fig == 'yes':
    plt.savefig('figure_s3.png', format='png', dpi=quality, bbox_inches = 'tight')
    
