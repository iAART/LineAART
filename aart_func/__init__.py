'''
Copyright (C) 2022, Alejandro Cardenas-Avendano, Alex Lupsasca & Hengrui Zhu
This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. 
If not, see https://www.gnu.org/licenses/.
'''
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
import h5py
import os
import argparse
import fileinput
import imageio.v2 as imageio

#For the analytical calculations
from numpy.lib.scimath import sqrt,  log
from numpy import arccos, arcsin,sin, cos
from scipy.special import ellipk, ellipeinc, ellipe
from scipy.special import ellipkinc as ellipf
from scipy.special import ellipj
from scipy.special import elliprj
from scipy.integrate import cumtrapz,quad

#Required for the wrapper for elliptic integral of the third kind
#import ctypes
#import numpy.ctypeslib as ctl

#For the lensing bands
from scipy.spatial import Delaunay
from matplotlib import path as paths

#Radon transformation
from skimage.transform import radon
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq, fftshift
from scipy import interpolate, optimize 

#Plotting 
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.spatial import delaunay_plot_2d
import matplotlib
cmap = plt.get_cmap('plasma') # This is the official colors of the BHs!

#Warnings flags
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
warnings.simplefilter("ignore", UserWarning)
np.seterr(divide='ignore', invalid='ignore')

#Various auxiliary functions
from aart_func.misc import *

#Elliptic function
from aart_func.ellipi_f import *

#Lensing Bands functions
import aart_func.lb_f as lb

def LensingBands(a, i):
    return lb.LensingBands(a, i)

#Raytracing
import aart_func.raytracing_f as rt

def RayTrace(aa,ii):
    return rt.RayTrace(aa,ii)

#Observer's Intensity
import aart_func.intensity_f as obsint

def Redshifts(aa,ii,SK,BBr,BBphi):
    return obsint.Redshifts(aa,ii,SK,BBr,BBphi)

def Intensities(aa,ii,SK,BBr,BBphi,Alpha_r,R_in0):
    return obsint.Intensities(aa,ii,SK,BBr,BBphi,Alpha_r,R_in0)

def RedshiftsAndIntensities(aa,ii,SK,BBr,BBphi,Alpha_r,R_in0):
    return obsint.RedshiftsAndIntensities(aa,ii,SK,BBr,BBphi,Alpha_r,R_in0)

#Line Profile 
import aart_func.lineprof_f as lp

def LineProfile(aa,ii,SK,BBr,BBphi,Alpha_r,R_in0,Es,Rmin,Rmax,Ebinwidth):
    return lp.LineProfile(aa,ii,SK,BBr,BBphi,Alpha_r,R_in0,Es,Rmin,Rmax,Ebinwidth)
