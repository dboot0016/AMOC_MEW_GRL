#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Geophysical Research Letters 06/2023
# Author script: A. Boot
# Contact person: Amber Boot (she/they; d.boot@uu.nl)

# Script for plotting Figure 2

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
vers = '1' # What Es-coupling equation to use: '1', '2' or '3'

save_fig = 'no' # Whether to save figure: 'yes' to save
quality = 300 # Quality figure

LW = 4 # Linewidth of line
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

#%% Load in branch 1
run1_1 = pd.read_csv(f'{Datadir}/'+'b.ref056_0_branch_1',sep='\s+',header=None)   
run2_1 = pd.read_csv(f'{Datadir}/'+'b.bio056_0_branch_1',sep='\s+',header=None)   
run3_1 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_0_branch_1',sep='\s+',header=None)   
run4_1 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_0_branch_1',sep='\s+',header=None)    

#%% Load in branch 2  
run2_2 = pd.read_csv(f'{Datadir}/'+'b.bio056_0_branch_2',sep='\s+',header=None)
run3_2 = pd.read_csv(f'{Datadir}/'+'b.es'+vers+'056_0_branch_2',sep='\s+',header=None)   
run4_2 = pd.read_csv(f'{Datadir}/'+'b.fca'+vers+'056_0_branch_2',sep='\s+',header=None)   

#%% Variables branch 1
# Ea
EA1_1 = run1_1[4]
EA2_1 = run2_1[4]
EA3_1 = run3_1[4]
EA4_1 = run4_1[4]

# Salinity box t
St1_1 = run1_1[6]*S_dim
St2_1 = run2_1[6]*S_dim
St3_1 = run3_1[6]*S_dim
St4_1 = run4_1[6]*S_dim

# Salinity box ts
Sts1_1 = run1_1[7]*S_dim
Sts2_1 = run2_1[7]*S_dim
Sts3_1 = run3_1[7]*S_dim
Sts4_1 = run4_1[7]*S_dim

# Salinity box n
Sn1_1 = run1_1[8]*S_dim
Sn2_1 = run2_1[8]*S_dim
Sn3_1 = run3_1[8]*S_dim
Sn4_1 = run4_1[8]*S_dim

# Salinity box s
Ss1_1 = run1_1[9]*S_dim
Ss2_1 = run2_1[9]*S_dim
Ss3_1 = run3_1[9]*S_dim
Ss4_1 = run4_1[9]*S_dim

# Pycnocline depth
D1_1 = run1_1[10]*D_dim
D2_1 = run2_1[10]*D_dim
D3_1 = run3_1[10]*D_dim
D4_1 = run4_1[10]*D_dim

# Atmospheric CO2 concentration 
CO21_1 = run1_1[23]*1e6
CO22_1 = run2_1[23]*1e6
CO23_1 = run3_1[23]*1e6
CO24_1 = run4_1[23]*1e6

#%% Variables branch 2 (similar as for branch 1)
EA2_2 = run2_2[4]
EA3_2 = run3_2[4]
EA4_2 = run4_2[4]

St2_2 = run2_2[6]*S_dim
St3_2 = run3_2[6]*S_dim
St4_2 = run4_2[6]*S_dim

Sts2_2 = run2_2[7]*S_dim
Sts3_2 = run3_2[7]*S_dim
Sts4_2 = run4_2[7]*S_dim

Sn2_2 = run2_2[8]*S_dim
Sn3_2 = run3_2[8]*S_dim
Sn4_2 = run4_2[8]*S_dim

Ss2_2 = run2_2[9]*S_dim
Ss3_2 = run3_2[9]*S_dim
Ss4_2 = run4_2[9]*S_dim

D2_2 = run2_2[10]*D_dim
D3_2 = run3_2[10]*D_dim
D4_2 = run4_2[10]*D_dim

CO22_2 = run2_2[23]*1e6
CO23_2 = run3_2[23]*1e6
CO24_2 = run4_2[23]*1e6

#%% Calculated variables branch 1
# Salinity box d
Sd1_1 = (S_dim*V0 - St1_1*(A*D1_1)-Sts1_1*(LxA*Ly*D1_1)-Ss1_1*Vs-Sn1_1*Vn)/(V0-Vn-Vs-LxA*Ly*D1_1-A*D1_1)
Sd2_1 = (S_dim*V0 - St2_1*(A*D2_1)-Sts2_1*(LxA*Ly*D2_1)-Ss2_1*Vs-Sn2_1*Vn)/(V0-Vn-Vs-LxA*Ly*D2_1-A*D2_1)
Sd3_1 = (S_dim*V0 - St3_1*(A*D3_1)-Sts3_1*(LxA*Ly*D3_1)-Ss3_1*Vs-Sn3_1*Vn)/(V0-Vn-Vs-LxA*Ly*D3_1-A*D3_1)
Sd4_1 = (S_dim*V0 - St4_1*(A*D4_1)-Sts4_1*(LxA*Ly*D4_1)-Ss4_1*Vs-Sn4_1*Vn)/(V0-Vn-Vs-LxA*Ly*D4_1-A*D4_1)

# AMOC strength
qn1_1 = 1e-6*(eta*alpha*(Tts-Tn)*(D1_1)**2+eta*beta*(Sn1_1-Sts1_1)*(D1_1)**2)
qn2_1 = 1e-6*(eta*alpha*(Tts-Tn)*(D2_1)**2+eta*beta*(Sn2_1-Sts2_1)*(D2_1)**2)
qn3_1 = 1e-6*(eta*alpha*(Tts-Tn)*(D3_1)**2+eta*beta*(Sn3_1-Sts3_1)*(D3_1)**2)
qn4_1 = 1e-6*(eta*alpha*(Tts-Tn)*(D4_1)**2+eta*beta*(Sn4_1-Sts4_1)*(D4_1)**2)

#%% Calculated variables branch 2
Sd3_2 = (S_dim*V0 - St3_2*(A*D3_2)-Sts3_2*(LxA*Ly*D3_2)-Ss3_2*Vs-Sn3_2*Vn)/(V0-Vn-Vs-LxA*Ly*D3_2-A*D3_2)
Sd4_2 = (S_dim*V0 - St4_2*(A*D4_2)-Sts4_2*(LxA*Ly*D4_2)-Ss4_2*Vs-Sn4_2*Vn)/(V0-Vn-Vs-LxA*Ly*D4_2-A*D4_2)
Sd2_2 = (S_dim*V0 - St2_2*(A*D2_2)-Sts2_2*(LxA*Ly*D2_2)-Ss2_2*Vs-Sn2_2*Vn)/(V0-Vn-Vs-LxA*Ly*D2_2-A*D2_2)

qn3_2 = 1e-6*(eta*alpha*(Tts-Tn)*(D3_2)**2+eta*beta*(Sn3_2-Sts3_2)*(D3_2)**2)
qn4_2 = 1e-6*(eta*alpha*(Tts-Tn)*(D4_2)**2+eta*beta*(Sn4_2-Sts4_2)*(D4_2)**2)
qn2_2 = 1e-6*(eta*alpha*(Tts-Tn)*(D2_2)**2+eta*beta*(Sn2_2-Sts2_2)*(D2_2)**2)

#%% Select data figure
# x-axis data (Ea)
x1 = EA1_1
x3 = EA2_1
x4 = EA2_2
x5 = EA3_1
x6 = EA3_2
x7 = EA4_1
x8 = EA4_2

# y-axis data
# CO2
if vary =='CO$_2$':
    y1 = CO21_1
    y3 = CO22_1
    y4 = CO22_2
    y5 = CO23_1
    y6 = CO23_2
    y7 = CO24_1
    y8 = CO24_2
    
    ymin = 0 # minimum y-axis
    ymax = 350 # maximum y-axis
    ymax2 = 350 # maximum y-axis
    unity = '[ppm]' # unit y-axis
    vary2 = vary # label y-axis

# AMOC strength    
elif vary =='D':
    y1 = (qn1_1*heaviside(qn1_1))
    y3 = (qn2_1*heaviside(qn2_1))
    y4 = qn2_2*heaviside(qn2_2)
    y5 = qn3_1*heaviside(qn3_1)
    y6 = qn3_2*heaviside(qn3_2)
    y7 = qn4_1*heaviside(qn4_1)
    y8 = qn4_2*heaviside(qn4_2)
    
    ymin = -0.75 # minimum y-axis
    ymax = 23 # maximum y-axis
    ymax2 = 23 # maximum y-axis
    unity = '[Sv]' # unit y-axis
    vary2 = 'AMOC' # label y-axis

#%%    
fig = plt.figure(figsize=(8, 6))

label1='REF'
label2='BIO'
label3='E$_s$ + BIO'
label4='E$_s$ + BIO + FCA'

if vers == '1':
    for i in np.arange(0,np.size(x1),1):
        if np.array(run1_1[1])[i]>0:
            LS='dotted' # Unstable branch
        else:
           LS='solid' # Stable branch   
        plt.plot(x1[i:i+2],y1[i:i+2],linewidth=LW,color='tab:blue',linestyle = LS,label='_nolegend_')
    
    plt.plot([-1e5,-2e5],[np.mean(y1),np.mean(y5)],linewidth=LW,color='tab:blue',label=label1)
    plt.plot([max(EA1_1),max(EA1_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:blue',zorder=1)
    
for i in np.arange(0,np.size(x3),1):
    if np.array(run2_1[1])[i]>0:
        LS='dotted' # Unstable branch
    else:
       LS='solid' # Stable branch   
    plt.plot(x3[i:i+2],y3[i:i+2],linewidth=LW,color='black',linestyle = LS,label='_nolegend_')

for i in np.arange(0,np.size(x4),1):
    if np.array(run2_2[1])[i]>0:
        LS='dotted' # Unstable branch
    else:
       LS='solid' # Stable branch    
    plt.plot(x4[i:i+2],y4[i:i+2],linewidth=LW,color='black',linestyle = LS,label='_nolegend_')

for i in np.arange(0,np.size(x5),1):
    if np.array(run3_1[1])[i]>0:
        LS='dotted' # Unstable branch
    else:
       LS='solid' # Stable branch    
    plt.plot(x5[i:i+2],y5[i:i+2],linewidth=LW,color='tab:orange',linestyle = LS,label='_nolegend_')

for i in np.arange(0,np.size(x6),1):
    if np.array(run3_2[1])[i]>0:
        LS='dotted' # Unstable branch
    else:
       LS='solid' # Stable branch    
    plt.plot(x6[i:i+2],y6[i:i+2],linewidth=LW,color='tab:orange',linestyle = LS,label='_nolegend_')

for i in np.arange(0,np.size(x7),1):
    if np.array(run4_1[1])[i]>0:
        LS='dotted' # Unstable branch
    else:
       LS='solid' # Stable branch   
    plt.plot(x7[i:i+2],y7[i:i+2],linewidth=LW,color='tab:green',linestyle = LS,label='_nolegend_')
    
for i in np.arange(0,np.size(x8),1):
    if np.array(run4_2[1])[i]>0:
        LS='dotted' # Unstable branch
    else:
       LS='solid'    # Stable branch
    plt.plot(x8[i:i+2],y8[i:i+2],linewidth=LW,color='tab:green',linestyle = LS,label='_nolegend_')

plt.plot([-1e5,-2e5],[np.mean(y1),np.mean(y5)],linewidth=LW,color='black',label=label2)
plt.plot([-1e5,-2e5],[np.mean(y1),np.mean(y5)],linewidth=LW,color='tab:orange',label=label3)
plt.plot([-1e5,-2e5],[np.mean(y1),np.mean(y5)],linewidth=LW,color='tab:green',label=label4)

# Saddle node locations
plt.plot([min(EA2_2),min(EA2_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='black',zorder=1)
plt.plot([max(EA2_1),max(EA2_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='black',zorder=1)
plt.plot([min(EA3_2),min(EA3_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:orange',zorder=1)
plt.plot([max(EA3_1),max(EA3_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:orange',zorder=1)
plt.plot([min(EA4_2),min(EA4_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:green',zorder=1)
plt.plot([max(EA4_1),max(EA4_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:green',zorder=1)

plt.xlabel(str(varx)+' ' +str(unitx),fontsize=FS-2)
plt.ylabel(str(vary2)+' ' +str(unity),fontsize=FS-2)

plt.xticks(fontsize=FS-4)
plt.yticks(fontsize=FS-4)    
plt.grid()

if vary=='D':
    plt.legend(fontsize=FS-6,loc=6)
plt.xlim([-0.22,0.5])
plt.ylim([ymin,ymax2])

if save_fig == 'yes':
    if vary=='D':
        plt.savefig('figure2_1_'+vers+'.png', format='png', dpi=quality, bbox_inches = 'tight')
    elif vary=='CO$_2$':
        plt.savefig('figure2_2_'+vers+'.png', format='png', dpi=quality, bbox_inches = 'tight')
