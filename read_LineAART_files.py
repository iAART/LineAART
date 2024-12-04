import numpy as np
import h5py

def rms(a):
    '''
    ISCO value
    :param a: BH spin
    '''
    Z1=1 + (1 - a**2)**(1/3) *((1 + a)**(1/3) + (1 - a)**(1/3))
    Z2=(3*a**2 + Z1**2)**(1/2)
    return (3 + Z2 - ((3 - Z1)*(3 + Z1 + 2*Z2))**(1/2))

def rh(a):
    '''
    "Outer horizon"
    :param a: BH spin
    '''
    return 1+np.sqrt(1-a**2)

def rph(a):
    '''
    PCO
    :param a: BH spin
    '''
    return 2*(1+np.cos((2/3)*np.arccos(-a)))

def namedrad(r,aa):
    if r == rms(aa):
        r='isco'
    elif r == rh(aa):
        r='horizon'
    elif r == rph(aa):
        r='pco'
    else:
        r=r
    return r

def rho_comp(rho):
    return 1-1/(rho+1) 

def rho_comp_inv(rhoC):
    return rhoC/(1-rhoC)

def SpecFldr(g0,g1,g2,aa,ii,path):
    return path+"SpectraData%s%s%s_a_%s_i_%s/"%(g0,g1,g2,aa,ii)

def LBFile(l0,l1,l2,aa,ii,path):
    return path+"LensingBands%s%s%s_a_%s_i_%s.h5"%(l0,l1,l2,aa,ii)

def RayFile(r0,r1,r2,aa,ii,path):
    return path+"Rays%s%s%s_a_%s_i_%s.h5"%(r0,r1,r2,aa,ii)

def gFile(g0,g1,g2,aa,ii,SK,BBr,BBphi,path):
    return SpecFldr(g0,g1,g2,aa,ii,path)+"gfactors%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s.h5" % (g0,g1,g2,aa,ii,SK,BBr,BBphi)

def IFile(g0,g1,g2,aa,ii,SK,BBr,BBphi,SigR,R0,path):
    return SpecFldr(g0,g1,g2,aa,ii,path)+"Intensities%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s_sig_%s_r0_%s.h5" % (g0,g1,g2,aa,ii,SK,BBr,BBphi,SigR,namedrad(R0,aa))

def LPFile(g0,g1,g2,aa,ii,SK,BBr,BBphi,SigR,R0,Es,Rmin,Rmax,path):
    return SpecFldr(g0,g1,g2,aa,ii,path)+"LineProfile%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s_sig_%s_r0_%s_Es_%s_rmin_%s_rmax_%s.h5" % (g0,g1,g2,aa,ii,SK,BBr,BBphi,SigR,namedrad(R0,aa),Es,namedrad(Rmin,aa),namedrad(Rmax,aa))