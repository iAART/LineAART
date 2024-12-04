from aart_func import *

inf=np.inf

def cbrt(x):
    '''
    Cubic root
    :param x: number to compute cbrt
    '''
    if x.imag==0:
        return np.cbrt(x)
    else:
        return x**(1/3)

def rms(a):
    '''
    ISCO value
    (Eq. B16 P1)
    :param a: BH spin
    '''
    Z1=1 + (1 - a**2)**(1/3) *((1 + a)**(1/3) + (1 - a)**(1/3))
    Z2=sqrt(3*a**2 + Z1**2)
    return (3 + Z2 - sqrt((3 - Z1)*(3 + Z1 + 2*Z2)))

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


def rho_comp(rho):
    return 1-1/(rho+1) 

def rho_comp_inv(rhoC):
    return rhoC/(1-rhoC)

def rho_by_rhoC(rhoC):    
    return (1/(1-rhoC)**2)


