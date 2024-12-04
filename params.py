from aart_func import *



########################################################
#params in this section: 
#  - are ignored when running Line-AART via in-line functions
#  - are used when running Line-AART via .py files in params_run folder
#  - can be change via params_run/Edit_params.py


#BH's Spin
spin_case=0.53
#Observer's inclination  
i_case=20

#Velocity Profile for the gas
#Sub-Kepleniarity param
sub_kep=0.9
#Radial velocity param
betar=0.9
#Angular velocity param
betaphi=0.9


#Disk Emissivity   I(r) = ( r^(-sigma) - r_0^(-sigma) ) Heaviside( r_0 - r )
#Emissivity fall off
sigma=3
#Emissivity cut off  
r_0=200

#Disk Emission Energy, in rest frame of gas 
E_s=7

#Disk Geometry, min and max disk radius
r_min=4.114523239768497
r_max=200


# energy bin width for spectrum
E_binwidth=0.05


########################################################
#params in this section: 
#  - are used when running Line-AART via in-line functions or .py files in params_run folder
#  - can be change via params_run/Edit_params.py

PrintStats=1

if PrintStats==True:
    print("\nThanks for using Line-AART")

#Limits for the image [M].
limits=15

#Resolution for the n=0 image [M], [radians]
drho0=0.001
dvarphi0=0.001

#Resolution for the n=1 image [M], [radians]
drho1=0.001
dvarphi1=0.001

#Resolution for the n=2 image [M], [radians]
drho2=0.001
dvarphi2=0.001

#Turn on/off the n=0,1,2 in LensingBands, RayTrace, and RedshiftFactors/Intensities
on_l0=1
on_l1=1
on_l2=1

on_r0=1
on_r1=1
on_r2=1

on_g0=1
on_g1=1
on_g2=1

N_on_l=[on_l0, on_l1, on_l2]
N_on_r=[on_r0, on_r1, on_r2]
N_on_g=[on_g0, on_g1, on_g2]


########################################################
#miscellaneous params

# Number of points in the critical curve   
npointsS=100 

#If equal to 1, the sizes of the grids will be equal and an image can be computed
#by summing the contributions    
p_image=0

#Current version is just implemented for equatorial models   
i_disk=90    
thetad=i_disk*np.pi/180


# Distance to the BH in meters (default: M87)
dBH=5.214795112e23  
# BH mass-to-distance ratio (default: 1/psi= 6.2e9 Kg)
psi=1.07473555940836 
#Observer's distance in units of M
D_obs=1e50   


Gc=6.67e-11 # G constant [m^3 kg^-1 s^-2]
cc= 2.99792458e8 # c constant [m/s]
Msc=1.988435e30 # Solar Mass [Kg]

MMkg= 6.2e9*psi*Msc # [Kg]
MM=MMkg *Gc/cc**2 # Mass of the BH in meters, i.e., for M87(psi*6.2*10^9) psi ("Best fit") Solar Masses 

# Size of the real image in meters
sizeim_Real=(limits)*MM 
#1 microarcsec in radians
muas_to_rad = np.pi/648000 *1e-6 
fov_Real=np.arctan(sizeim_Real/(dBH))/muas_to_rad #muas
#print("FOV= ",np.round(2*fov,2),"muas")


###############################################################################
#results path in this section: 
#  - can be change via params_run/Edit_params.py

#Path where the results will be stored
path="./Results_Examples_Part2/"

# Create a directory for the results
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
    print("A Results directory  was created to store the results")

    
###############################################################################    

'''
MIT license
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.
'''



