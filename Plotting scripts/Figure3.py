#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Geophysical Research Letters 06/2023
# Author script: A. Boot
# Contact person: Amber Boot (she/they; d.boot@uu.nl)

# Script for plotting Figure 3

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

#%% Set up plotting domain
y1 = 0.995 + 0.001
y2 = 0.995 + 0.003
y3 = 0.995 + 0.005
y4 = 0.995 + 0.007
y5 = 0.995 + 0.009

y6 = 0.995 + 0.011
y7 = 0.995 + 0.013
y8 = 0.995 + 0.015
y9 = 0.995 + 0.017
y10 = 0.995 + 0.019

y11 = 0.995 + 0.021
y12 = 0.995 + 0.023
y13 = 0.995 + 0.025
y14 = 0.995 + 0.027
y15 = 0.995 + 0.029

#%%
run15_1 = pd.read_csv(f'{Datadir}/'+'b.bio056_m2000_branch_1',sep='\s+',header=None)   
run10_1 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_m2000_branch_1',sep='\s+',header=None)   
run5_1 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_m2000_branch_1',sep='\s+',header=None)    

run15_2 = pd.read_csv(f'{Datadir}/'+'b.bio056_m2000_branch_2',sep='\s+',header=None)   
run10_2 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_m2000_branch_2',sep='\s+',header=None)   
run5_2 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_m2000_branch_2',sep='\s+',header=None)   

#%% Load in branch 1
run14_1 = pd.read_csv(f'{Datadir}/'+'b.bio056_0_branch_1',sep='\s+',header=None)   
run9_1 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_0_branch_1',sep='\s+',header=None)   
run4_1 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_0_branch_1',sep='\s+',header=None)    
 
run4_2 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_0_branch_2',sep='\s+',header=None)
run14_2 = pd.read_csv(f'{Datadir}/'+'b.bio056_0_branch_2',sep='\s+',header=None)
run9_2 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_0_branch_2',sep='\s+',header=None)

#%% Load in branch 1
run12_1 = pd.read_csv(f'{Datadir}/'+'b.bio056_4000_branch_1',sep='\s+',header=None)   
run7_1 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_4000_branch_1',sep='\s+',header=None)    
run2_1 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_4000_branch_1',sep='\s+',header=None)    
  
run2_2 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_4000_branch_2',sep='\s+',header=None)
run12_2 = pd.read_csv(f'{Datadir}/'+'b.bio056_4000_branch_2',sep='\s+',header=None)
run7_2 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_4000_branch_2',sep='\s+',header=None)

#%% Load in branch 1
run11_1 = pd.read_csv(f'{Datadir}/'+'b.bio056_8000_branch_1',sep='\s+',header=None)   
run6_1 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_8000_branch_1',sep='\s+',header=None)    
run1_1 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_8000_branch_1',sep='\s+',header=None)    
  
run1_2 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_8000_branch_2',sep='\s+',header=None)
run11_2 = pd.read_csv(f'{Datadir}/'+'b.bio056_8000_branch_2',sep='\s+',header=None)
run6_2 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_8000_branch_2',sep='\s+',header=None)

#%% Load in branch 1
run13_1 = pd.read_csv(f'{Datadir}/'+'b.bio056_2000_branch_1',sep='\s+',header=None)   
run8_1 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_2000_branch_1',sep='\s+',header=None)    
run3_1 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_2000_branch_1',sep='\s+',header=None)    
 
run13_2 = pd.read_csv(f'{Datadir}/'+'b.bio056_2000_branch_2',sep='\s+',header=None)
run8_2 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_2000_branch_2',sep='\s+',header=None)
run3_2 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_2000_branch_2',sep='\s+',header=None)

#%% Select Ea data
EA1_1 = run1_1[4]
EA2_1 = run2_1[4]
EA3_1 = run3_1[4]
EA4_1 = run4_1[4]
EA5_1 = run5_1[4]

EA6_1 = run6_1[4]
EA7_1 = run7_1[4]
EA8_1 = run8_1[4]
EA9_1 = run9_1[4]
EA10_1 = run10_1[4]

EA11_1 = run11_1[4]
EA12_1 = run12_1[4]
EA13_1 = run13_1[4]
EA14_1 = run14_1[4]
EA15_1 = run15_1[4]

EA1_2 = run1_2[4]
EA2_2 = run2_2[4]
EA3_2 = run3_2[4]
EA4_2 = run4_2[4]
EA5_2 = run5_2[4]

EA6_2 = run6_2[4]
EA7_2 = run7_2[4]
EA8_2 = run8_2[4]
EA9_2 = run9_2[4]
EA10_2 = run10_2[4]

EA11_2 = run11_2[4]
EA12_2 = run12_2[4]
EA13_2 = run13_2[4]
EA14_2 = run14_2[4]
EA15_2 = run15_2[4]

#%% Select atmospheric pCO2
CO21_1 = run1_1[23]*1e6
CO22_1 = run2_1[23]*1e6
CO23_1 = run3_1[23]*1e6
CO24_1 = run4_1[23]*1e6
CO25_1 = run5_1[23]*1e6

CO26_1 = run6_1[23]*1e6
CO27_1 = run7_1[23]*1e6
CO28_1 = run8_1[23]*1e6
CO29_1 = run9_1[23]*1e6
CO210_1 = run10_1[23]*1e6

CO211_1 = run11_1[23]*1e6
CO212_1 = run12_1[23]*1e6
CO213_1 = run13_1[23]*1e6
CO214_1 = run14_1[23]*1e6
CO215_1 = run15_1[23]*1e6

CO21_2 = run1_2[23]*1e6
CO22_2 = run2_2[23]*1e6
CO23_2 = run3_2[23]*1e6
CO24_2 = run4_2[23]*1e6
CO25_2 = run5_2[23]*1e6

CO26_2 = run6_2[23]*1e6
CO27_2 = run7_2[23]*1e6
CO28_2 = run8_2[23]*1e6
CO29_2 = run9_2[23]*1e6
CO210_2 = run10_2[23]*1e6

CO211_2 = run11_2[23]*1e6
CO212_2 = run12_2[23]*1e6
CO213_2 = run13_2[23]*1e6
CO214_2 = run14_2[23]*1e6
CO215_2 = run15_2[23]*1e6

#%% Determine Ea value of saddle nodes
# ES + BIO + FCA
# On branch
x1_1 = max(EA1_1)
x2_1 = max(EA2_1)
x3_1 = max(EA3_1)
x4_1 = max(EA4_1)
x5_1 = max(EA5_1)

# Off branch
x1_2 = min(EA1_2)
x2_2 = min(EA2_2)
x3_2 = min(EA3_2)
x4_2 = min(EA4_2)
x5_2 = min(EA5_2)

# ES + BIO
# On branch
x6_1 = max(EA6_1)
x7_1 = max(EA7_1)
x8_1 = max(EA8_1)
x9_1 = max(EA9_1)
x10_1 = max(EA10_1)

# Off branch
x6_2 = min(EA6_2)
x7_2 = min(EA7_2)
x8_2 = min(EA8_2)
x9_2 = min(EA9_2)
x10_2 = min(EA10_2)

# BIO
# On branch
x11_1 = max(EA11_1)
x12_1 = max(EA12_1)
x13_1 = max(EA13_1)
x14_1 = max(EA14_1)
x15_1 = max(EA15_1)

# Off branch
x11_2 = min(EA11_2)
x12_2 = min(EA12_2)
x13_2 = min(EA13_2)
x14_2 = min(EA14_2)
x15_2 = min(EA15_2)

#%% Figure 3a and d
Y=np.zeros(100)+1
xline=np.linspace(-1.5,1.5,np.size(Y))

X_out = np.linspace(10,11,100)
Y_out = np.linspace(10,11,100)

plot_lines = []
plot_lines2 = []

fig, ax = plt.subplots(figsize=(8,6))

L5,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:purple',label = '-2000 PgC')
L4,=plt.plot(X_out,Y_out,linewidth=LW,color='black',label = '0 PgC')
L3,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:green',label = '+2000 PgC')
L2,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:orange',label = '+4000 PgC')
L1,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:blue',label = '+8000 PgC')

l1=plt.scatter(X_out,Y_out,s=MS/2,color='black',label = 'On')
l2=plt.scatter(X_out,Y_out,s=MS/2,marker='s',color='black',label = 'Off')

plot_lines.append([l1, l2])
plot_lines2.append([L5, L4, L3, L2, L1])

plt.scatter(x1_1,y1,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x2_1,y2,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x3_1,y3,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x4_1,y4,s=MS,color='black',label='_nolegend_')
plt.scatter(x5_1,y5,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x1_2,y1,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x2_2,y2,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x3_2,y3,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x4_2,y4,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x5_2,y5,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x1_1,x1_2],[y1,y1],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x2_1,x2_2],[y2,y2],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x3_1,x3_2],[y3,y3],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x4_1,x4_2],[y4,y4],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x5_1,x5_2],[y5,y5],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.scatter(x6_1,y6,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x7_1,y7,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x8_1,y8,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x9_1,y9,s=MS,color='black',label='_nolegend_')
plt.scatter(x10_1,y10,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x6_2,y6,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x7_2,y7,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x8_2,y8,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x9_2,y9,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x10_2,y10,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x6_1,x6_2],[y6,y6],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x7_1,x7_2],[y7,y7],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x8_1,x8_2],[y8,y8],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x9_1,x9_2],[y9,y9],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x10_1,x10_2],[y10,y10],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.scatter(x11_1,y11,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x12_1,y12,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x13_1,y13,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x14_1,y14,s=MS,color='black',label='_nolegend_')
plt.scatter(x15_1,y15,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x11_2,y11,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x12_2,y12,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x13_2,y13,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x14_2,y14,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x15_2,y15,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x11_1,x11_2],[y11,y11],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x12_1,x12_2],[y12,y12],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x13_1,x13_2],[y13,y13],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x14_1,x14_2],[y14,y14],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x15_1,x15_2],[y15,y15],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.plot(xline,Y+0.005,'k',linewidth=0.75)
plt.plot(xline,Y+0.015,'k',linewidth=0.75)
plt.plot(xline,Y+0.025,'k',linewidth=0.75)

#plt.title('Bistability window',fontsize=FS)
plt.xticks(fontsize=FS-4)
plt.xlim([-0.25, 0.7])
plt.ylim([0.995, 1.025])
plt.yticks([1,1.01,1.02],['FCA','E$_s$','BIO'],fontsize=FS-2)
plt.xlabel('E$_a$ [Sv]',fontsize=FS-2)
ax.xaxis.grid(True, which="minor", ls="--")
ax.xaxis.grid(True, which="major", ls="-")

if save_fig == 'on':
    plt.savefig('figure_3_1_' + vers + '.png', format='png', dpi=quality, bbox_inches = 'tight')
    

#%%
plot_lines = []
plot_lines2 = []

#%% Determine CO2 values corresponding to saddle nodes
# ES + BIO + FCA
# On branch
x1_1 = CO21_1[np.where(max(EA1_1)==EA1_1)[0][0]]
x2_1 = CO22_1[np.where(max(EA2_1)==EA2_1)[0][0]]
x3_1 = CO23_1[np.where(max(EA3_1)==EA3_1)[0][0]]
x4_1 = CO24_1[np.where(max(EA4_1)==EA4_1)[0][0]]
x5_1 = CO25_1[np.where(max(EA5_1)==EA5_1)[0][0]]

# Off branch
x1_2 = CO21_2[np.where(min(EA1_2)==EA1_2)[0][0]]
x2_2 = CO22_2[np.where(min(EA2_2)==EA2_2)[0][0]]
x3_2 = CO23_2[np.where(min(EA3_2)==EA3_2)[0][0]]
x4_2 = CO24_2[np.where(min(EA4_2)==EA4_2)[0][0]]
x5_2 = CO25_2[np.where(min(EA5_2)==EA5_2)[0][0]]

# ES + BIO
# On branch
x6_1 = CO26_1[np.where(max(EA6_1)==EA6_1)[0][0]]
x7_1 = CO27_1[np.where(max(EA7_1)==EA7_1)[0][0]]
x8_1 = CO28_1[np.where(max(EA8_1)==EA8_1)[0][0]]
x9_1 = CO29_1[np.where(max(EA9_1)==EA9_1)[0][0]]
x10_1 = CO210_1[np.where(max(EA10_1)==EA10_1)[0][0]]

# Off branch
x6_2 = CO26_2[np.where(min(EA6_2)==EA6_2)[0][0]]
x7_2 = CO27_2[np.where(min(EA7_2)==EA7_2)[0][0]]
x8_2 = CO28_2[np.where(min(EA8_2)==EA8_2)[0][0]]
x9_2 = CO29_2[np.where(min(EA9_2)==EA9_2)[0][0]]
x10_2 = CO210_2[np.where(min(EA10_2)==EA10_2)[0][0]]

# BIO
# On branch
x11_1 = CO211_1[np.where(max(EA11_1)==EA11_1)[0][0]]
x12_1 = CO212_1[np.where(max(EA12_1)==EA12_1)[0][0]]
x13_1 = CO213_1[np.where(max(EA13_1)==EA13_1)[0][0]]
x14_1 = CO214_1[np.where(max(EA14_1)==EA14_1)[0][0]]
x15_1 = CO215_1[np.where(max(EA15_1)==EA15_1)[0][0]]

# Off branch
x11_2 = CO211_2[np.where(min(EA11_2)==EA11_2)[0][0]]
x12_2 = CO212_2[np.where(min(EA12_2)==EA12_2)[0][0]]
x13_2 = CO213_2[np.where(min(EA13_2)==EA13_2)[0][0]]
x14_2 = CO214_2[np.where(min(EA14_2)==EA14_2)[0][0]]
x15_2 = CO215_2[np.where(min(EA15_2)==EA15_2)[0][0]]

#%% Figure 3b, e
Y=np.zeros(100)+1
xline=np.linspace(0,1550,np.size(Y))

fig, ax = plt.subplots(figsize=(8,6))

L5,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:purple',label = '-2000 PgC')
L4,=plt.plot(X_out,Y_out,linewidth=LW,color='black',label = '0 PgC')
L3,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:green',label = '+2000 PgC')
L2,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:orange',label = '+4000 PgC')
L1,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:blue',label = '+8000 PgC')

l1=plt.scatter(X_out,Y_out,s=MS/2,color='black',label = 'On')
l2=plt.scatter(X_out,Y_out,s=MS/2,marker='s',color='black',label = 'Off')

plot_lines.append([l1, l2])
plot_lines2.append([L5, L4, L3, L2, L1])

plt.scatter(x1_1,y1,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x2_1,y2,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x3_1,y3,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x4_1,y4,s=MS,color='black',label='_nolegend_')
plt.scatter(x5_1,y5,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x1_2,y1,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x2_2,y2,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x3_2,y3,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x4_2,y4,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x5_2,y5,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x1_1,x1_2],[y1,y1],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x2_1,x2_2],[y2,y2],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x3_1,x3_2],[y3,y3],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x4_1,x4_2],[y4,y4],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x5_1,x5_2],[y5,y5],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.scatter(x6_1,y6,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x7_1,y7,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x8_1,y8,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x9_1,y9,s=MS,color='black',label='_nolegend_')
plt.scatter(x10_1,y10,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x6_2,y6,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x7_2,y7,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x8_2,y8,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x9_2,y9,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x10_2,y10,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x6_1,x6_2],[y6,y6],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x7_1,x7_2],[y7,y7],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x8_1,x8_2],[y8,y8],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x9_1,x9_2],[y9,y9],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x10_1,x10_2],[y10,y10],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.scatter(x11_1,y11,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x12_1,y12,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x13_1,y13,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x14_1,y14,s=MS,color='black',label='_nolegend_')
plt.scatter(x15_1,y15,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x11_2,y11,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x12_2,y12,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x13_2,y13,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x14_2,y14,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x15_2,y15,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x11_1,x11_2],[y11,y11],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x12_1,x12_2],[y12,y12],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x13_1,x13_2],[y13,y13],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x14_1,x14_2],[y14,y14],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x15_1,x15_2],[y15,y15],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.plot(xline,Y+0.005,'k',linewidth=0.75)
plt.plot(xline,Y+0.015,'k',linewidth=0.75)
plt.plot(xline,Y+0.025,'k',linewidth=0.75)

#plt.title('Bistability window',fontsize=FS)
ax.set_xscale('log')
plt.xticks(fontsize=FS-4)
plt.xlim([0, 1500])
plt.ylim([0.995, 1.025])
plt.yticks([1,1.01,1.02],['FCA','E$_s$','BIO'],fontsize=FS-2)
plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
ax.xaxis.grid(True, which="minor", ls="--")
ax.xaxis.grid(True, which="major", ls="-")


if save_fig == 'on':
    plt.savefig('figure_3_2_' + vers + '.png', format='png', dpi=quality, bbox_inches = 'tight')
    
#%% Es - CO2 coupling to select Es value corresponding to saddle ndoes
if vers == '2':
    ln_a = 0.099
    ln_b = 0.069
    ln_c = 0
    ln_d = 1
elif vers == '3':
    ln_a = 0.1
    ln_b = 0.75
    ln_c = 250
    ln_d = 300

ES1_1 = ln_a+(ln_b*np.log((x1_1+ln_c)/ln_d))
ES2_1 = ln_a+(ln_b*np.log((x2_1+ln_c)/ln_d))
ES3_1 = ln_a+(ln_b*np.log((x3_1+ln_c)/ln_d))
ES4_1 = ln_a+(ln_b*np.log((x4_1+ln_c)/ln_d))
ES5_1 = ln_a+(ln_b*np.log((x5_1+ln_c)/ln_d))

ES6_1 = ln_a+(ln_b*np.log((x6_1+ln_c)/ln_d))
ES7_1 = ln_a+(ln_b*np.log((x7_1+ln_c)/ln_d))
ES8_1 = ln_a+(ln_b*np.log((x8_1+ln_c)/ln_d))
ES9_1 = ln_a+(ln_b*np.log((x9_1+ln_c)/ln_d))
ES10_1 = ln_a+(ln_b*np.log((x10_1+ln_c)/ln_d))

ES11_1 = ln_a+(ln_b*np.log((x11_1+ln_c)/ln_d))
ES12_1 = ln_a+(ln_b*np.log((x12_1+ln_c)/ln_d))
ES13_1 = ln_a+(ln_b*np.log((x13_1+ln_c)/ln_d))
ES14_1 = ln_a+(ln_b*np.log((x14_1+ln_c)/ln_d))
ES15_1 = ln_a+(ln_b*np.log((x15_1+ln_c)/ln_d))

ES1_2 = ln_a+(ln_b*np.log((x1_2+ln_c)/ln_d))
ES2_2 = ln_a+(ln_b*np.log((x2_2+ln_c)/ln_d))
ES3_2 = ln_a+(ln_b*np.log((x3_2+ln_c)/ln_d))
ES4_2 = ln_a+(ln_b*np.log((x4_2+ln_c)/ln_d))
ES5_2 = ln_a+(ln_b*np.log((x5_2+ln_c)/ln_d))

ES6_2 = ln_a+(ln_b*np.log((x6_2+ln_c)/ln_d))
ES7_2 = ln_a+(ln_b*np.log((x7_2+ln_c)/ln_d))
ES8_2 = ln_a+(ln_b*np.log((x8_2+ln_c)/ln_d))
ES9_2 = ln_a+(ln_b*np.log((x9_2+ln_c)/ln_d))
ES10_2 = ln_a+(ln_b*np.log((x10_2+ln_c)/ln_d))

ES11_2 = ln_a+(ln_b*np.log((x11_2+ln_c)/ln_d))
ES12_2 = ln_a+(ln_b*np.log((x12_2+ln_c)/ln_d))
ES13_2 = ln_a+(ln_b*np.log((x13_2+ln_c)/ln_d))
ES14_2 = ln_a+(ln_b*np.log((x14_2+ln_c)/ln_d))
ES15_2 = ln_a+(ln_b*np.log((x15_2+ln_c)/ln_d))

#%% Select the data
# ES + FCA + BIO
# On branch
x1_1 = ES1_1
x2_1 = ES2_1
x3_1 = ES3_1
x4_1 = ES4_1
x5_1 = ES5_1

# Off branch
x1_2 = ES1_2
x2_2 = ES2_2
x3_2 = ES3_2
x4_2 = ES4_2
x5_2 = ES5_2

# ES + BIO
# On branch
x6_1 = ES6_1
x7_1 = ES7_1
x8_1 = ES8_1
x9_1 = ES9_1
x10_1 = ES10_1

# Off branch
x6_2 = ES6_2
x7_2 = ES7_2
x8_2 = ES8_2
x9_2 = ES9_2
x10_2 = ES10_2

# BIO
# On branch
x11_1 = ES11_1
x12_1 = ES12_1
x13_1 = ES13_1
x14_1 = ES14_1
x15_1 = ES15_1

# Off branch
x11_2 = ES11_2
x12_2 = ES12_2
x13_2 = ES13_2
x14_2 = ES14_2
x15_2 = ES15_2

#%% Figure 3c, f
Y=np.zeros(100)+1
xline=np.linspace(-100,1550,np.size(Y))

fig, ax = plt.subplots(figsize=(8,6))

L5,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:purple',label = '-2000 PgC')
L4,=plt.plot(X_out,Y_out,linewidth=LW,color='black',label = '0 PgC')
L3,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:green',label = '+2000 PgC')
L2,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:orange',label = '+4000 PgC')
L1,=plt.plot(X_out,Y_out,linewidth=LW,color='tab:blue',label = '+8000 PgC')

l1=plt.scatter(X_out,Y_out,s=MS/2,color='black',label = 'On')
l2=plt.scatter(X_out,Y_out,s=MS/2,marker='s',color='black',label = 'Off')

plot_lines.append([l1, l2])
plot_lines2.append([L5, L4, L3, L2, L1])

plt.scatter(x1_1,y1,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x2_1,y2,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x3_1,y3,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x4_1,y4,s=MS,color='black',label='_nolegend_')
plt.scatter(x5_1,y5,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x1_2,y1,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x2_2,y2,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x3_2,y3,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x4_2,y4,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x5_2,y5,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x1_1,x1_2],[y1,y1],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x2_1,x2_2],[y2,y2],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x3_1,x3_2],[y3,y3],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x4_1,x4_2],[y4,y4],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x5_1,x5_2],[y5,y5],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.scatter(x6_1,y6,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(x7_1,y7,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(x8_1,y8,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(x9_1,y9,s=MS,color='black',label='_nolegend_')
plt.scatter(x10_1,y10,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(x6_2,y6,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(x7_2,y7,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(x8_2,y8,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(x9_2,y9,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(x10_2,y10,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot([x6_1,x6_2],[y6,y6],linewidth=LW,linestyle='dashed',color='tab:blue',label='_nolegend_')
plt.plot([x7_1,x7_2],[y7,y7],linewidth=LW,linestyle='dashed',color='tab:orange',label='_nolegend_')
plt.plot([x8_1,x8_2],[y8,y8],linewidth=LW,linestyle='dashed',color='tab:green',label='_nolegend_')
plt.plot([x9_1,x9_2],[y9,y9],linewidth=LW,linestyle='dashed',color='black',label='_nolegend_')
plt.plot([x10_1,x10_2],[y10,y10],linewidth=LW,linestyle='dashed',color='tab:purple',label='_nolegend_')

plt.scatter(0.56,y11,s=MS,color='tab:blue',label='_nolegend_')
plt.scatter(0.56,y12,s=MS,color='tab:orange',label='_nolegend_')
plt.scatter(0.56,y13,s=MS,color='tab:green',label='_nolegend_')
plt.scatter(0.56,y14,s=MS,color='black',label='_nolegend_')
plt.scatter(0.56,y15,s=MS,color='tab:purple',label='_nolegend_')

plt.scatter(0.56,y11,s=MS,marker='s',color='tab:blue',label='_nolegend_')
plt.scatter(0.56,y12,s=MS,marker='s',color='tab:orange',label='_nolegend_')
plt.scatter(0.56,y13,s=MS,marker='s',color='tab:green',label='_nolegend_')
plt.scatter(0.56,y14,s=MS,marker='s',color='black',label='_nolegend_')
plt.scatter(0.56,y15,s=MS,marker='s',color='tab:purple',label='_nolegend_')

plt.plot(xline,Y+0.005,'k',linewidth=0.75)
plt.plot(xline,Y+0.015,'k',linewidth=0.75)
plt.plot(xline,Y+0.025,'k',linewidth=0.75)

plt.xticks(fontsize=FS-4)

if vers == '2':
    plt.xlim([0.15, 0.7])
elif vers == '3':
    plt.xlim([-0.1,1.2])
    
plt.ylim([0.995, 1.025])
plt.yticks([1,1.01,1.02],['FCA','E$_s$','BIO'],fontsize=FS-2)
plt.xlabel('E$_s$ [Sv]',fontsize=FS-2)
ax.xaxis.grid(True, which="minor", ls="--")
ax.xaxis.grid(True, which="major", ls="-")

legend1 = plt.legend(handles=plot_lines[0],fontsize=FS-7,loc='upper right')
plt.gca().add_artist(legend1)

legend2 = plt.legend(handles=plot_lines2[0],fontsize=FS-7,loc='upper left')
plt.gca().add_artist(legend2)

if save_fig == 'on':
    plt.savefig('figure_3_3_' + vers + '.png', format='png', dpi=quality, bbox_inches = 'tight')

