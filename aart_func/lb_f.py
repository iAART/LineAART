from aart_func import *
from params import * 

#The inverse of the radial potential Eq. 5 P1. 
def OverRint(r,a,lamb,eta):
    return 1/sqrt((r**2 + a**2 - a*lamb)**2 - (r**2 - 2*r + a**2)*(eta + (lamb - a)**2))
    
def ApparentBH(s,a,thetao,alpha,beta,m,sign,distance=D_obs):

    rho=np.sqrt(alpha**2+beta**2)
    varphi=np.arctan2(beta,alpha)

    alpha= s*rho*np.cos(varphi)
    beta= s*rho*np.sin(varphi)
    
    # Photon conserved quantities
    # Eqs. (55 P1)
    lamb = -alpha*np.sin(thetao) 
    eta = (alpha**2 - a**2)*np.cos(thetao)**2+beta**2
    
    DeltaTheta = 1/2 *(1 - (eta + lamb**2)/a**2) # Eq. (11 P1)
    
    # Roots of angular potentail
    uP = DeltaTheta + sqrt(DeltaTheta**2 + eta/a**2) # Eq. (11 P1)
    uM = DeltaTheta - sqrt(DeltaTheta**2 + eta/a**2) # Eq. (11 P1)
    
    #Outer horizon of the BH
    rP = 1 + np.sqrt(1 - a**2)

    #We perform the integration Eq. 8a using Eq. 20 P1 in order to get zero. 
    return sqrt(-uM*a**2)*quad(OverRint, rP,distance,args=(a,lamb,eta),points=[rP,10000])[0] +sign*ellipf(np.arcsin(np.cos(thetao)/np.sqrt(uP)), uP/uM) -2*m*ellipk(uP/uM) 

def ApparentBH0(s,a,alpha,beta,distance=D_obs):

    rho=np.sqrt(alpha**2+beta**2)
    varphi=np.arctan2(beta,alpha)

    alpha= s*rho*np.cos(varphi)
    beta= s*rho*np.sin(varphi)
    
    # Photon conserved quantities
    # Eqs. (55 P1)
    eta = (alpha**2 +beta**2 - a**2)
    
    #Outer horizon of the BH
    rP = 1 + np.sqrt(1 - a**2)
    
    return sqrt(eta)*quad(OverRint, rP,distance,args=(a,0,eta),points=[rP,1000])[0]-ellipk(-a**2/eta)

# def OverRint(rinv,a,lamb,eta):
#     r=1/rinv
#     return (-1/rinv**2)*1/sqrt((r**2 + a**2 - a*lamb)**2 - (r**2 - 2*r + a**2)*(eta + (lamb - a)**2))

# def ApparentBH(s,a,thetao,alpha,beta,m,sign,distance=D_obs):

#     rho=np.sqrt(alpha**2+beta**2)
#     varphi=np.arctan2(beta,alpha)

#     alpha= s*rho*np.cos(varphi)
#     beta= s*rho*np.sin(varphi)
    
#     # Photon conserved quantities
#     # Eqs. (55 P1)
#     lamb = -alpha*np.sin(thetao) 
#     eta = (alpha**2 - a**2)*np.cos(thetao)**2+beta**2
    
#     DeltaTheta = 1/2 *(1 - (eta + lamb**2)/a**2) # Eq. (11 P1)
    
#     # Roots of angular potentail
#     uP = DeltaTheta + sqrt(DeltaTheta**2 + eta/a**2) # Eq. (11 P1)
#     uM = DeltaTheta - sqrt(DeltaTheta**2 + eta/a**2) # Eq. (11 P1)
    
#     #Outer horizon of the BH
#     rP = 1 + np.sqrt(1 - a**2)
    
#     #distinace_inverse=1/D_obs
#     distinace_inverse=0

#     #We perform the integration Eq. 8a using Eq. 20 P1 in order to get zero. 
#     return sqrt(-uM*a**2)*quad(OverRint, 1/rP,distinace_inverse,args=(a,lamb,eta))[0] +sign*ellipf(np.arcsin(np.cos(thetao)/np.sqrt(uP)), uP/uM) -2*m*ellipk(uP/uM) 

def nlayers(s,a,thetao,thetad,alpha,betai,mbar):
    '''
    Computes the boundary of the nth lensing band
    :param s: Type of ray. s>1 (Rays arrive outside) and s<1 (Rays arrive inside)
              the critical curve. We have to solve for this parameter. 
    :param a: BH's spin (-1<0<1)
    :param thetao: Observers angle [degrees]
    :param thetad: Angle of the disk [degrees]
    :param alpha: Bardeen's coordinate alpha
    :param betai: Bardeen's coordinate beta
    :param mbar: Label of the observed rings 

    :return: 
    '''
    alpha= s*alpha
    beta= s*betai

    #Angular Turning points encountered along the trajectory
    m= mbar + np.heaviside(betai,0) # Eq. (82 P1)

    # Photon conserved quantities
    # Eqs. (55 P1)
    lamb = -alpha*np.sin(thetao) 
    eta = (alpha**2 - a**2)*np.cos(thetao)**2+beta**2
    
    nutheta=np.sign(betai)
    
    # Radial Roots and Integrals 
    AAA = a**2 - eta - lamb**2 # Eq. (79 P3)
    BBB = 2*(eta + (lamb - a)**2) # Eq. (80 P3)
    CCC = -a**2 *eta # Eq. (81 P3)
    
    P = -(AAA**2/12) - CCC # Eq. (85 P3)
    Q = -(AAA/3) *((AAA/6)**2 - CCC) - BBB**2/8  # Eq. (86 P3)

    Delta3 = -4 *P**3 - 27*Q**2 # Eq. (92 P3)
    
    xi0 = np.real(cbrt(-(Q/2) + sqrt(-(Delta3/108))) + cbrt(-(Q/2) - sqrt(-(Delta3/108))) - AAA/3) # Eq. (87 P3)
    z = sqrt(xi0/2) # Eq. (94 P3)
   
    r1 = -z - sqrt(-(AAA/2) - z**2 + BBB/(4*z)) # Eq. (95a P3)
    r2 = -z + sqrt(-(AAA/2) - z**2 + BBB/(4*z)) # Eq. (95b P3)
    r3 = z - sqrt(-(AAA/2) - z**2 - BBB/(4*z))  # Eq. (95c P3)
    r4 = z + sqrt(-(AAA/2) - z**2 - BBB/(4*z))  # Eq. (95d P3)
    
    DeltaTheta = 1/2 *(1 - (eta + lamb**2)/a**2) # Eq. (19 P3)
    
    # Roots of angular potentail
    # Eqs. (B9 P3)
    uP = DeltaTheta + sqrt(DeltaTheta**2 + eta/a**2) # Eq. (19 P3)
    uM = DeltaTheta - sqrt(DeltaTheta**2 + eta/a**2) # Eq. (19 P3)
    
    # Eqs. (B9 P3)
    r21 = r2 - r1 
    r31 = r3 - r1
    r32 = r3 - r2
    r41 = r4 - r1
    r42 = r4 - r2
    r43 = r4 - r3
    
    # Outer and inner horizons
    # Eqs. (2 P3)
    rP = 1 + np.sqrt(1 - a**2)
    rM = 1 - np.sqrt(1 - a**2)
    
    # Eqs. (B10 P3)
    a1=sqrt(-(r43**2/4))
    # Eqs. (B10 P3)
    b1=(r3 + r4)/2
    
    #Elliptic Parameter
    # Eqs. (B13 P3)
    k = (r32*r41)/(r31*r42)

    AA = np.real(sqrt(a1**2 + (b1 - r2)**2)) # Eqs. (B56 P3)
    BB = np.real(sqrt(a1**2 + (b1 - r1)**2)) # Eqs. (B56 P3)

    # This parameter is real and less the unity
    k3 = np.real(((AA + BB)**2 - r21**2)/(4*AA*BB)) # Eqs. (B59 P3)
    
    # Eqs. (20 P2)
    Gtheta = 1/(sqrt(-uM)*a)*(2*m*ellipk(uP/uM) -nutheta*ellipf(np.arcsin(np.cos(thetao)/np.sqrt(uP)), uP/uM) + nutheta*(-1)**m*ellipf(np.arcsin(np.cos(thetad)/np.sqrt(uP)), uP/uM))
    
    if s>1:
        # Eqs.  (A10 P2)
        Q1= 4/sqrt(r31*r42)*ellipf(np.arcsin(sqrt(r31/r41)), k)
        return Q1-Gtheta
    else:
        if k3<1:
            # Eqs. (A11 P2)
            Q2=1/sqrt(AA*BB)*(ellipf(np.arccos((AA - BB)/(AA + BB)), k3) - ellipf(np.arccos((AA *(rP - r1) - BB*(rP - r2))/(AA*(rP - r1) + BB*(rP - r2))), k3))
        else:
            Q2=np.nan
    
        return Q2-Gtheta


def nlayers0(s,a,alpha,betai,mbar):
    '''
    Computes the boundary of the nth lensing band
    :param s: Type of ray. s>1 (Rays arrive outside) and s<1 (Rays arrive inside)
              the critical curve. We have to solve for this parameter. 
    :param a: BH's spin (-1<0<1)
    :param thetao: Observers angle [degrees]
    :param thetad: Angle of the disk [degrees]
    :param alpha: Bardeen's coordinate alpha
    :param betai: Bardeen's coordinate beta
    :param mbar: Label of the observed rings 

    :return: 
    '''
    alpha= s*alpha
    beta= s*betai

    # Photon conserved quantities
    # Eqs. (55 P1)
    eta = alpha**2 - a**2+beta**2
        
    # Radial Roots and Integrals 
    AAA = a**2 - eta # Eq. (79 P3)
    BBB = 2*(eta + a**2) # Eq. (80 P3)
    CCC = -a**2 *eta # Eq. (81 P3)
    
    P = -(AAA**2/12) - CCC # Eq. (85 P3)
    Q = -(AAA/3) *((AAA/6)**2 - CCC) - BBB**2/8  # Eq. (86 P3)

    Delta3 = -4 *P**3 - 27*Q**2 # Eq. (92 P3)
    
    xi0 = np.real(cbrt(-(Q/2) + sqrt(-(Delta3/108))) + cbrt(-(Q/2) - sqrt(-(Delta3/108))) - AAA/3) # Eq. (87 P3)
    z = sqrt(xi0/2) # Eq. (94 P3)
   
    r1 = -z - sqrt(-(AAA/2) - z**2 + BBB/(4*z)) # Eq. (95a P3)
    r2 = -z + sqrt(-(AAA/2) - z**2 + BBB/(4*z)) # Eq. (95b P3)
    r3 = z - sqrt(-(AAA/2) - z**2 - BBB/(4*z))  # Eq. (95c P3)
    r4 = z + sqrt(-(AAA/2) - z**2 - BBB/(4*z))  # Eq. (95d P3)
    
    # Eqs. (B9 P3)
    r21 = r2 - r1 
    r31 = r3 - r1
    r32 = r3 - r2
    r41 = r4 - r1
    r42 = r4 - r2
    r43 = r4 - r3
    
    # Outer and inner horizons
    # Eqs. (2 P3)
    rP = 1 + np.sqrt(1 - a**2)
    rM = 1 - np.sqrt(1 - a**2)
    
    # Eqs. (B10 P3)
    a1=sqrt(-(r43**2/4))
    # Eqs. (B10 P3)
    b1=(r3 + r4)/2
    
    #Elliptic Parameter
    # Eqs. (B13 P3)
    k = (r32*r41)/(r31*r42)

    AA = np.real(sqrt(a1**2 + (b1 - r2)**2)) # Eqs. (B56 P3)
    BB = np.real(sqrt(a1**2 + (b1 - r1)**2)) # Eqs. (B56 P3)

    # This parameter is real and less the unity
    k3 = np.real(((AA + BB)**2 - r21**2)/(4*AA*BB)) # Eqs. (B59 P3)
    
    #
    Gtheta = ((1+2*mbar)*ellipk(-a**2/eta)/sqrt(eta)).real
      
    if s>1:
        # Eqs.  (A10 P2)
        Q1= 4/sqrt(r31*r42)*ellipf(np.arcsin(sqrt(r31/r41)), k)
        return Q1-Gtheta
    else:
        if k3<1:
            # Eqs. (A11 P2)
            Q2=1/sqrt(AA*BB)*(ellipf(np.arccos((AA - BB)/(AA + BB)), k3) - ellipf(np.arccos((AA *(rP - r1) - BB*(rP - r2))/(AA*(rP - r1) + BB*(rP - r2))), k3))
        else:
            Q2=np.nan
    
        return Q2-Gtheta        
    
    
def spacedmarks(x, y, Nmarks):

    """
    Computes the arch-length
    :param x: x point
    :param y: y point
    :param Nmarks: Number of marks 

    :returns: position of the x and y markers
    """
    dydx = np.gradient(y, x[0],edge_order=2)
    dxdx = np.gradient(x, x[0],edge_order=2)
    arclength = cumtrapz(sqrt(dydx**2 + dxdx**2), initial=0)
    marks = np.linspace(0, max(arclength), Nmarks)
    markx = np.interp(marks, arclength, x)
    marky = np.interp(markx, x, y)
    return markx, marky
    

def CritCurve(a,angle):
    """
    Computes the critical curve (aka the black hole shadow)
    :param a: spin of black hole
    :param angle: angle of the observer in degrees

    :returns: contour of the critical curve on the observer plane
    """
    thetao = angle * np.pi/180
    
    #Critical roots of the radial potential (Eq. 5 P2) given by (Eq. 40 P2)
    rM = 2*(1 + np.cos(2/3 *np.arccos(-(a))))
    rP = 2*(1 + np.cos(2/3 *np.arccos(a)))
    
    #Define the region where \tilde{r} lies
    r=np.linspace(rM,rP,int(1e7))
    
    #Eqs. 38 & 39 P2
    lam = a + r/a *(r - (2 *(r**2 - 2*r + a**2))/(r - 1))
    eta = r**3/a**2 *((4*(r**2 - 2*r + a**2))/(r - 1)**2 - r)

    #Eqs. 55 P2
    alpha=-lam/np.sin(thetao)
    # This is beta squared actually, but we need to check where it is positive
    beta=eta + a**2 *np.cos(thetao)**2 - lam**2*np.tan(thetao)**(-2)
    
    mask=np.where(beta>0)
    r=r[mask]
    
    rmin=min(r)+1e-12
    rmax=max(r)-1e-12
    
    r=np.linspace(rmin,rmax,int(1e6))
        
    #Eqs. 38 & 39 P2
    lam = a + r/a *(r - (2 *(r**2 - 2*r + a**2))/(r - 1))
    eta = r**3/a**2 *((4*(r**2 - 2*r + a**2))/(r - 1)**2 - r)

    #Eqs. 55 P2
    alpha=-lam/np.sin(thetao)
    beta=eta + a**2 *np.cos(thetao)**2 - lam**2*np.tan(thetao)**(-2)

    return alpha, sqrt(beta)

def CritCurve0(a):
    """
    Computes the critical curve (aka the black hole shadow)
    :param a: spin of black hole
    :param angle: angle of the observer in degrees

    :returns: contour of the critical curve on the observer plane
    """
      
    xx=(9 - 9 * a ** 2 + sqrt(3) * a * sqrt(-27 + 18 * a ** 2 + a ** 4))**(1/3)
    r=(1-(a**2-3)/(3 ** (1/3) * xx)+xx/(3**(2/3))).real
    eta = r**3/a**2 *((4*(r**2 - 2*r + a**2))/(r - 1)**2 - r).real
    b=sqrt(eta + a**2).real 
    
    ang=np.linspace(1e-12,np.pi-1e-12,npointsS)
    
    alpha=b*np.cos(ang)
    beta=b*np.sin(ang)
    
    return alpha, beta

def round_up_to_even(f):
    return int(np.ceil(f / 2.) * 2)

def in_hull(p, hull):
    """
    Check if points in p are inside the concave hull
    """
    concave=paths.Path(hull)
    return concave.contains_points(p)
    

def grid_mask_polar(hull, hull2, drho, dvarphi, limits, force_lims = False):
    """
    Create polar grid on the observer plane
    :param hull: marks inner edge of the lensing band
    :param hull2: marks outer edge of the lensing band
    :param drho: spacing of the radial direction on the screen
    :param dvarphi: spacing of the angular direction on the screen
    :param limits: specify maximum screen radius
    :param force_lims: specify a limit or find the optimal limit for the grid

    :returns: a polar grid in Cartesian coordinates and other parameters
    """
    # We set the limits of rho list
    if force_lims == False:
        
        rho_min = np.sqrt(min(hull2[:,0]**2 + hull2[:,1]**2)) - 5 * dvarphi
        rho_max = np.sqrt(max(hull2[:,0]**2 + hull2[:,1]**2)) + 5 * dvarphi
        
        lims = np.ceil(max(np.abs(rho_min), np.abs(rho_max)))
        
        if (lims >= limits):
            lims = limits
        
    else:
        lims = limits
        
    # We find the rho and varphi list we want keeping spacing less than drho and dvarphi respectively
    Nrho=int(round_up_to_even(1/drho))
    rho=rho_comp(lims)*np.linspace(1/(Nrho), 1, Nrho)
        
    Nvarphi=int(round_up_to_even(2*np.pi/dvarphi))
    varphi=2*np.pi*(np.linspace(1/(Nvarphi), 1, Nvarphi)-1/(4*(Nvarphi)))
   
    # We create grid of Cartesian points from rho and varphi lists   
    mesh = np.array(np.meshgrid(rho, varphi))
    grid = mesh.T.reshape(-1, 2)
    
    rho, varphi = zip(*grid)
    rho=tuple(rho_comp_inv(np.array(rho)))
    grid_new = list(zip(rho * np.cos(varphi), rho * np.sin(varphi)))

    
    # We check the points belong to the lensing bands
    mask1 = np.invert(in_hull(grid_new, hull))
    mask2 = in_hull(grid_new, hull2)

    indexes = mask2 * mask1

    return np.array(grid_new), Nrho, Nvarphi, indexes, lims


def hullsN(a,thetao, alpha, beta, smin=0.5, smax=100,limi0=0.99,lime0=1.01,limi1=0.999,lime1=1.001,limi2=0.9999,lime2=1.001):
    """
    :param a: black hole spin
    :thetao: observer inclination in radians
    :param alpha: x coordinate on screen
    :param beta: y corrdinate on screen
    
    :returns: hulls the inner and outer edges of the lensing bands
    """
    
    data = (np.append(alpha, alpha[::-1]), np.append(beta, -beta[::-1]))
    
    N0_on=bool(N_on_l[0])
    N1_on=bool(N_on_l[1])
    N2_on=bool(N_on_l[2])
    
    if N0_on == True: 
        points_0i=np.zeros([data[0].size, 2])
        # Do this clockwise, to maximize the number of points contained in the hull. 
        points_0e=np.array([[-limits,-limits],[limits,-limits],[limits,limits],[-limits, limits]])
    else:
        points_0i=np.array([[0,0]])
        points_0e=np.array([[0,0]])
  
    if N1_on == True: 
        points_1i=np.zeros([data[0].size, 2])
        points_1e=np.zeros([data[0].size, 2])
    else:
        points_1i=np.array([[0,0]])
        points_1e=np.array([[0,0]])
        
    if N2_on == True: 
        points_2i=np.zeros([data[0].size, 2])
        points_2e=np.zeros([data[0].size, 2])
    else:
        points_2i=np.array([[0,0]])
        points_2e=np.array([[0,0]])


    for i in range(data[0].size):
        
        if N0_on == True:
            
            if thetao==0:
                m1=optimize.root(ApparentBH0,limi0, args=(a,data[0][i],data[1][i],D_obs))
            else:
                if data[1][i]>=0:
                    m1=optimize.root(ApparentBH, limi0, args=(a,thetao,data[0][i],data[1][i],1,1,D_obs))
                else:
                    m1=optimize.root(ApparentBH, limi0, args=(a,thetao,data[0][i],data[1][i],0,-1,D_obs))

            points_0i[i]=m1.x[0]*np.array([data[0][i],data[1][i]])
            
            
        if N1_on == True:

            m1=optimize.root(nlayers, limi1, args=(a,thetao,thetad,data[0][i],data[1][i],1))
            m2=optimize.root(nlayers, lime1, args=(a,thetao,thetad,data[0][i],data[1][i],1))

            points_1i[i]=limi1*m1.x[0]*np.array([data[0][i],data[1][i]])
            points_1e[i]=lime1*m2.x[0]*np.array([data[0][i],data[1][i]])

            if m1.x[0]<smin:
                points_1i[i]=limi1*smin*np.array([data[0][i],data[1][i]])
            else:
                points_1i[i]=limi1*m1.x[0]*np.array([data[0][i],data[1][i]]) 

            if m2.x[0]>smax:
                points_1e[i]=lime1*smax*np.array([data[0][i],data[1][i]])
            else:
                points_1e[i]=lime1*m2.x[0]*np.array([data[0][i],data[1][i]])
                
        if N2_on == True:

            m1=optimize.root(nlayers, limi2, args=(a,thetao,thetad,data[0][i],data[1][i],2))
            m2=optimize.root(nlayers, lime2, args=(a,thetao,thetad,data[0][i],data[1][i],2))

            points_2i[i]=limi2*m1.x[0]*np.array([data[0][i],data[1][i]])
            points_2e[i]=lime2*m2.x[0]*np.array([data[0][i],data[1][i]])

    return points_0i, points_0e, points_1i, points_1e, points_2i, points_2e

def LensingBands(a, i):
    """
    :param a: black hole spin
    :thetao: observer inclination in radians
    
    :return: LensingBand file containing the following for each lensing band n that is turned on
    :output hull_ni: inner edges of lensing bands in cartesian coordinates
    :output hull_ne: outer edges of lensing bands in cartesian coordinates
    :output gridn: screen coordinates (rho,varphi) of points 
    :output maskn: mask for points inside lensing band
    :output limn: maximum screen radius rho
    :output Nrhon: number of screen radii in gridn 
    :output Nrhon:number of screen angles gridn Nphin
    
    """
    
    if PrintStats==True:
        print('Computing Lensing Bands')
    # "Computing the grids for each lensing bands"
    if i==0:
        alpha_critc, beta_critc = CritCurve0(a)
        #critc = CritCurve0(a)
        #alpha_critc, beta_critc = spacedmarks(critc[0], critc[1], npointsS)
    else:     
        critc = CritCurve(a, i)
        alpha_critc, beta_critc = spacedmarks(critc[0], critc[1], npointsS)
        #alpha_critc, beta_critc = CritCurve(a, i)
    
    thetao=i*np.pi/180

    hull_0i, hull_0e, hull_1i, hull_1e, hull_2i, hull_2e = hullsN(a,thetao,alpha_critc, beta_critc)
    
    N0_on=bool(N_on_l[0])
    N1_on=bool(N_on_l[1])
    N2_on=bool(N_on_l[2])

    if (p_image == 1):
        
        if N0_on == True:
            supergrid0, Nrho0, Nvarphi0, mask0, lim0 = grid_mask_polar(hull_0i,hull_0e,drho0,dvarphi0,limits,True)
            if PrintStats==True:
                print("Number of points in the n=0 grid", supergrid0.shape[0])
            
        if N1_on == True:
            supergrid1, Nrho1, Nvarphi1, mask1, lim1 = grid_mask_polar(hull_1i,hull_1e,drho1,dvarphi1,limits,True)
            if PrintStats==True:
                print("Number of points in the n=1 grid", supergrid1.shape[0])
    
        if N2_on == True:
            supergrid2, Nrho2, Nvarphi2, mask2, lim2 = grid_mask_polar(hull_2i,hull_2e,drho2,dvarphi2,limits,True)
            if PrintStats==True:
                print("Number of points in the n=2 grid", supergrid2.shape[0])

        
    else:
        
        if N0_on == True:
            supergrid0, Nrho0, Nvarphi0, mask0, lim0 = grid_mask_polar(hull_0i,hull_0e,drho0,dvarphi0,limits)
            if PrintStats==True:
                print("Number of points in the n=0 grid", supergrid0.shape[0])
            
        if N1_on == True:
            supergrid1, Nrho1, Nvarphi1, mask1, lim1 = grid_mask_polar(hull_1i,hull_1e,drho1,dvarphi1,limits)
            if PrintStats==True:
                print("Number of points in the n=1 grid", supergrid1.shape[0])
     
        if N2_on == True:
            supergrid2, Nrho2, Nvarphi2, mask2, lim2 = grid_mask_polar(hull_2i,hull_2e,drho2,dvarphi2,limits)
            if PrintStats==True:
                print("Number of points in the n=2 grid", supergrid2.shape[0])


    filename=path+"LensingBands%s%s%s_a_%s_i_%s.h5"%(int(N0_on),int(N1_on),int(N2_on),a,i)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('alpha', data=alpha_critc)
    h5f.create_dataset('beta', data=beta_critc)
    
    if N0_on == True:
        h5f.create_dataset('hull_0i', data=hull_0i)
        h5f.create_dataset('hull_0e', data=hull_0e)
        h5f.create_dataset('grid0', data=supergrid0)
        h5f.create_dataset('Nrho0', data=np.array([Nrho0]))
        h5f.create_dataset('Nvarphi0', data=np.array([Nvarphi0]))
        h5f.create_dataset('mask0', data=mask0)
        h5f.create_dataset('lim0', data=np.array([lim0]))
        
    if N1_on == True:
        h5f.create_dataset('hull_1i', data=hull_1i)
        h5f.create_dataset('hull_1e', data=hull_1e)
        h5f.create_dataset('grid1', data=supergrid1)
        h5f.create_dataset('Nrho1', data=np.array([Nrho1]))
        h5f.create_dataset('Nvarphi1', data=np.array([Nvarphi1]))
        h5f.create_dataset('mask1', data=mask1)
        h5f.create_dataset('lim1', data=np.array([lim1]))
        
    if N2_on == True:
        h5f.create_dataset('hull_2i', data=hull_2i)
        h5f.create_dataset('hull_2e', data=hull_2e)
        h5f.create_dataset('grid2', data=supergrid2)
        h5f.create_dataset('Nrho2', data=np.array([Nrho2]))
        h5f.create_dataset('Nvarphi2', data=np.array([Nvarphi2]))
        h5f.create_dataset('mask2', data=mask2)
        h5f.create_dataset('lim2', data=np.array([lim2]))
    
    
    h5f.create_dataset('N_on_l', data=np.array([N_on_l[0],N_on_l[1],N_on_l[2]]))
    h5f.close()

    if PrintStats==True:
        print("File ",filename," created.")