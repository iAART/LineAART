from aart_func import *
from params import * 

def Delta(r,a):
    """
    Calculates the Kerr metric function \Delta(t)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return r**2-2*r+a**2

def PIF(r,a):
    """
    Calculates PI(r) (Eq. B6 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (r**2+a**2)**2-a**2*Delta(r,a)

def urbar(r,a):
    """
    Calculates the r (contravariant) component of the four velocity for radial infall
    (Eq. B34b P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return -np.sqrt(2*r*(r**2+a**2))/(r**2)

def Omegabar(r,a):
    """
    Calculates the angular velocity of the radial infall
    (Eq. B32a P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (2*a*r)/PIF(r,a)

def Omegahat(r,a,laux):
    """
    Calculates the angular velocity of the sub-Keplerian orbit
    (Eq. B39 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (a+(1-2/r)*(laux-a))/(PIF(r,a)/(r**2)-(2*a*laux)/r)

def uttilde(r, a,urT,OT):
    """
    Calculates the t (contravariant) component of the general four velocity
    (Eq. B52 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param urT: r (contravariant) component of the general four velocity
    :param OT: Angular velocity of the general four velocity
    """
    return np.sqrt((1 + urT**2*r**2/Delta(r,a))/(1-(r**2+a**2)*OT**2-(2/r)*(1-a*OT)**2))

def Ehat(r,a,laux):
    """
    Calculates the orbital energy of the sub-Keplerian flow
    (Eq. B44a P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param laux: sub-Keplerian specific angular momentum
    """
    return np.sqrt(Delta(r,a)/(PIF(r,a)/(r**2)-(4*a*laux)/r-(1-2/r)*laux**2))

def nuhat(r,a,laux,Ehataux):
    """
    Calculates the radial velocity of the sub-Keplerian flow
    (Eq. B45 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param laux: sub-Keplerian specific angular momentum
    :param Ehataux: sub-Keplerian orbital energy
    """
    return r/Delta(r,a)*np.sqrt(np.abs(PIF(r,a)/(r**2)-(4*a*laux)/r-(1-2/r)*laux**2-Delta(r,a)/(Ehataux**2)))

def lhat(r,a,Sub_Kep):
    """
    Calculates the rspecific angular momentum of the sub-Keplerian flow
    (Eq. B44b P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param Sub_Kep: sub-Keplerian factor
    """
    return Sub_Kep*(r**2+a**2-2*a*np.sqrt(r))/(np.sqrt(r)*(r-2)+a)

def Rint(r,a,lamb,eta):
    """
    Evaluates the "radial potential", for calculating the redshift factor for infalling material
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: carter constant

    :return: radial potential evaluated at the source
    """
    #Eqns (P2 5)
    return (r**2 + a**2 - a*lamb)**2 - (r**2 - 2*r + a**2)*(eta + (lamb - a)**2)

def gDisk(r,a,b,lamb,eta,Sub_Kep,BetaR,BetaPhi):
    """
    Calculates the redshift factor for a photon outside the inner-most stable circular orbit(isco) (assume circular orbit)
    (Eq. B13 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: Carter constant
    :param Sub_Kep: sub-Keplerian factor
    :param BetaR, BetaPhi: raidal and angular disk motion mixing parameters

    :return: the redshift factor associated with the ray
    """

    OH=Omegahat(r,a,lhat(r,a,Sub_Kep))
    OT=OH+(1-BetaPhi)*(Omegabar(r,a)-OH)
    ur=(1-BetaR)*urbar(r,a)
    ut=uttilde(r,a,ur,OT)
    uphi=ut*OT
    
    return 1/(ut*(1-b*np.sign(ur)*sqrt(np.abs(Rint(r,a,lamb,eta)*ur**2))/Delta(r,a)/ut-lamb*uphi/ut))

def gGas(r,a,b,lamb,eta,Sub_Kep,BetaR,BetaPhi):
    """
    Calculates the redshift factor for a photon inside the isco (assume infalling orbit)
    (Eq. B13 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param b: sign for the redshift
    :param lamb: angular momentum
    :param eta: carter constant
    :params Sub_Kep,BetaR, BetaPhi: disk motion parameters

    :return: the redshift factor associated with the ray
    """
    #Calculate radius of the inner-most stable circular orbit
    isco=rms(a)

    lms=lhat(isco,a,Sub_Kep)
    OH=Omegahat(r,a,lms)
    OT=OH+(1-BetaPhi)*(Omegabar(r,a)-OH)

    Ems=Ehat(isco,a,lms)
    urhat=-Delta(r,a)/(r**2)*nuhat(r, a, lms ,Ems)*Ems
    ur=urhat+(1-BetaR)*(urbar(r,a)-urhat)
    ut=uttilde(r,a,ur,OT)
    uphi=OT*ut

    return 1/(ut*(1-b*np.sign(ur)*sqrt(np.abs(Rint(r,a,lamb,eta)*ur**2))/Delta(r,a)/ut-lamb*uphi/ut))


def gfactorf(a,i,grid,mask,redshift_sign,rs,Sub_Kep,BetaR,BetaPhi):
    """
    Calculate the redshift factor
    :param a: black hole spin
    :param i: observer inclination
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param rs: source radius
    :params Sub_Kep,BetaR, BetaPhi: disk motion parameters

    :return: redshift factor at each point.

    """
    
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    
    thetao=i*np.pi/180
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    gfact = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]
   

    isco=rms(a)
    r_p = rh(a)
    
    gfact[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco],Sub_Kep,BetaR,BetaPhi)
    gfact[(r_p<rs)&(rs<isco)]= gGas(rs[(r_p<rs)&(rs<isco)],a,redshift_sign[(r_p<rs)&(rs<isco)],lamb[(r_p<rs)&(rs<isco)],eta[(r_p<rs)&(rs<isco)],Sub_Kep,BetaR,BetaPhi)
    gfact[rs<=r_p] = np.nan
    
    gs = np.nan*np.zeros(mask.shape)
    gs[mask] = gfact
    return(gs)


def intensity(rs,g,alpha_r,r_int0):
    """
    Calculate the observed intensity of each pixel   
    :param rs: rs source radius
    :param g: redshift factor
    :alpha_r: power of source intensity function
    """
    if r_int0==inf:
        In=(rs**-alpha_r)* g**3
    else:
        In=(rs**-alpha_r - r_int0**-alpha_r) * np.heaviside(r_int0, rs) * g**3
    return(In)



def Redshifts(a,i,Sub_Kep,BetaR,BetaPhi):
    """
    Calculate the redshift factor from LensingBands and Rays files
    :param a: black hole spin
    :param i: observer inclination
    :params Sub_Kep, BetaR, BetaPhi: disk motion parameters

    :return: redshift factor at each point on on bserver screen.

    """
    
    if PrintStats==True:
        print("Computing Redshift")


    if N_on_l[0]<N_on_g[0] or N_on_l[1]< N_on_g[1] and N_on_l[2]< N_on_g[2]:

        if PrintStats==True:
            print("Missing lensing band files")
            print("N_on_l ",N_on_l)
            print("N_on_r ",N_on_g)

    else:
        if N_on_r[0]<N_on_g[0] or N_on_r[1]< N_on_g[1] and N_on_r[2]< N_on_g[2]:

            if PrintStats==True:
                print("Missing ray-tracing files")
                print("N_on_r ",N_on_l)
                print("N_on_g ",N_on_g)

        else:
            

            fn=path+"LensingBands%s%s%s_a_%s_i_%s.h5"%(N_on_l[0],N_on_l[1],N_on_l[2],a,i)

            N0_on=bool(N_on_g[0])
            N1_on=bool(N_on_g[1])
            N2_on=bool(N_on_g[2])

            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                supergrid0=h5f['grid0'][:]
                mask0=h5f['mask0'][:]
            else:
                supergrid0=np.array([])
                mask0=[]      

            if N1_on == True:
                supergrid1=h5f['grid1'][:]
                mask1=h5f['mask1'][:]
            else:
                supergrid1=np.array([])
                mask1=[]             

            if N2_on == True:    
                supergrid2=h5f['grid2'][:]
                mask2=h5f['mask2'][:]
            else:
                supergrid2=np.array([])
                mask2=[]     

            h5f.close()


            fn=path+"Rays%s%s%s_a_%s_i_%s.h5"%(N_on_r[0],N_on_r[1],N_on_r[2],a,i)


            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                rs0=h5f['rs0'][:]
                sign0=h5f['sign0'][:]

            if N1_on == True:
                rs1=h5f['rs1'][:]
                sign1=h5f['sign1'][:]

            if N2_on == True:    
                rs2=h5f['rs2'][:]
                sign2=h5f['sign2'][:]

            h5f.close()

            if N0_on == True:
                if PrintStats==True:
                    print("Calculating redshift factors for n=0,  # of points", supergrid0.shape[0])
                i_g0 = gfactorf(a,i,supergrid0,mask0,sign0,rs0,Sub_Kep,BetaR,BetaPhi) 

            if N1_on == True:
                if PrintStats==True:
                    print("Calculating redshift factors for n=1,  # of points", supergrid1.shape[0])
                i_g1 = gfactorf(a,i,supergrid1,mask1,sign1,rs1,Sub_Kep,BetaR,BetaPhi)

            if N2_on == True:
                if PrintStats==True:
                    print("Calculating redshift factors for n=2,  # of points", supergrid2.shape[0])
                i_g2 = gfactorf(a,i,supergrid2,mask2,sign2,rs2,Sub_Kep,BetaR,BetaPhi)

            gf_fldr="SpectraData%s%s%s_a_%s_i_%s/"%(N_on_l[0],N_on_l[1],N_on_l[2],a,i)
            os.makedirs(path+gf_fldr, exist_ok=True)
                
            filename=path+gf_fldr+"gfactors%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s.h5"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi)

            h5f = h5py.File(filename, 'w')

            if N0_on == True:
                h5f.create_dataset('gf0', data=i_g0)

            if N1_on == True:
                h5f.create_dataset('gf1', data=i_g1)

            if N1_on == True:
                h5f.create_dataset('gf2', data=i_g2)

            h5f.close()
            
            if PrintStats==True:
                print("File ",filename," created.")
            
            

            
def Intensities(a,i,Sub_Kep,BetaR,BetaPhi,Alpha_r,R_int0):
    """
    Calculate the intensities factor from LensingBands, Rays, and gFactors files
    :param a: black hole spin
    :param i: observer inclination
    :params Sub_Kep, BetaR, BetaPhi: disk motion parameters

    :return: intensity at each point on on bserver screen.

    """
    
    if PrintStats==True:
        print("Computing Intensity")


    if N_on_l[0]<N_on_g[0] or N_on_l[1]< N_on_g[1] and N_on_l[2]< N_on_g[2]:

        if PrintStats==True:
            print("Missing lensing band files")
            print("N_on_l ",N_on_l)
            print("N_on_r ",N_on_g)

    else:
        if N_on_r[0]<N_on_g[0] or N_on_r[1]< N_on_g[1] and N_on_r[2]< N_on_g[2]:

            if PrintStats==True:
                print("Missing ray-tracing files")
                print("N_on_r ",N_on_l)
                print("N_on_g ",N_on_g)

        else:

            fn=path+"LensingBands%s%s%s_a_%s_i_%s.h5"%(N_on_l[0],N_on_l[1],N_on_l[2],a,i)

            N0_on=bool(N_on_g[0])
            N1_on=bool(N_on_g[1])
            N2_on=bool(N_on_g[2])

            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                supergrid0=h5f['grid0'][:]
                mask0=h5f['mask0'][:]     

            if N1_on == True:
                supergrid1=h5f['grid1'][:]
                mask1=h5f['mask1'][:]

            if N2_on == True:    
                supergrid2=h5f['grid2'][:]
                mask2=h5f['mask2'][:]

            h5f.close()


            fn=path+"Rays%s%s%s_a_%s_i_%s.h5"%(N_on_r[0],N_on_r[1],N_on_r[2],a,i)


            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                rs0=h5f['rs0'][:]

            if N1_on == True:
                rs1=h5f['rs1'][:]

            if N2_on == True:    
                rs2=h5f['rs2'][:]

            h5f.close()
           
            gf_fldr="SpectraData%s%s%s_a_%s_i_%s/"%(N_on_l[0],N_on_l[1],N_on_l[2],a,i)
            
            fn=path+ gf_fldr + "gfactors%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s.h5"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi)


            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                gf0=h5f['gf0'][:]

            if N1_on == True:
                gf1=h5f['gf1'][:]

            if N2_on == True:    
                gf2=h5f['gf2'][:]

            h5f.close() 
            
   
            #r_p = rh(a)

            if N0_on == True:
                if PrintStats==True:
                    print("Calculating intensity factors for n=0,  # of points", supergrid0.shape[0])
                I_0 = np.nan*np.zeros(rs0.shape[0])
                I_0[mask0] = intensity(rs0[mask0],gf0[mask0],Alpha_r,R_int0) 

            if N1_on == True:
                if PrintStats==True:
                    print("Calculating intensity factors for n=1,  # of points", supergrid1.shape[0])
                I_1 = np.nan*np.zeros(rs1.shape[0])
                I_1[mask1] = intensity(rs1[mask1],gf1[mask1],Alpha_r,R_int0)

            if N2_on == True:
                if PrintStats==True:
                    print("Calculating intensity factors for n=2,  # of points", supergrid2.shape[0])
                I_2 = np.nan*np.zeros(rs2.shape[0])
                I_2[mask2] = intensity(rs2[mask2],gf2[mask2],Alpha_r,R_int0)
        
            if R_int0 == rms(a):
                Rint02='isco'
                
            elif R_int0 == rh(a):
                Rint02='horizon'
                
            elif R_int0 == rph(a):
                Rint02='pco'    
                
            else:
                Rint02=R_int0
            
            filename=path+gf_fldr+"Intensities%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s_sig_%s_r0_%s.h5" %(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi,Alpha_r,Rint02)

            h5f = h5py.File(filename, 'w')

            if N0_on == True:
                h5f.create_dataset('Io0', data=I_0)

            if N1_on == True:
                h5f.create_dataset('Io1', data=I_1)

            if N1_on == True:
                h5f.create_dataset('Io2', data=I_2)

            h5f.close()
            
            if PrintStats==True:
                print("File ",filename," created.")   
            

def RedshiftsAndIntensities(a,i,Sub_Kep,BetaR,BetaPhi,Alpha_r,R_int0):
    """
    Calculate the redshift factors and intensities factor from LensingBands, Rays, and gFactors files
    :param a: black hole spin
    :param i: observer inclination
    :params Sub_Kep, BetaR, BetaPhi: disk motion parameters

    :return: redshift factor and intensity at each point on on bserver screen.

    """
    
    if PrintStats==True:
        print("Computing Redshift and Intensity")


    if N_on_l[0]<N_on_g[0] or N_on_l[1]< N_on_g[1] and N_on_l[2]< N_on_g[2]:

        if PrintStats==True:
            print("Missing lensing band files")
            print("N_on_l ",N_on_l)
            print("N_on_r ",N_on_g)

    else:
        if N_on_r[0]<N_on_g[0] or N_on_r[1]< N_on_g[1] and N_on_r[2]< N_on_g[2]:

            if PrintStats==True:
                print("Missing ray-tracing files")
                print("N_on_r ",N_on_l)
                print("N_on_g ",N_on_g)

        else:

            fn=path+"LensingBands%s%s%s_a_%s_i_%s.h5"%(N_on_l[0],N_on_l[1],N_on_l[2],a,i)

            N0_on=bool(N_on_g[0])
            N1_on=bool(N_on_g[1])
            N2_on=bool(N_on_g[2])

            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                supergrid0=h5f['grid0'][:]
                mask0=h5f['mask0'][:]   

            if N1_on == True:
                supergrid1=h5f['grid1'][:]
                mask1=h5f['mask1'][:]       

            if N2_on == True:    
                supergrid2=h5f['grid2'][:]
                mask2=h5f['mask2'][:]    

            h5f.close()


            fn=path+"Rays%s%s%s_a_%s_i_%s.h5"%(N_on_r[0],N_on_r[1],N_on_r[2],a,i)


            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                rs0=h5f['rs0'][:]
                sign0=h5f['sign0'][:]

            if N1_on == True:
                rs1=h5f['rs1'][:]
                sign1=h5f['sign1'][:]

            if N2_on == True:    
                rs2=h5f['rs2'][:]
                sign2=h5f['sign2'][:]

            h5f.close()
            
            
            r_p = rh(a)
            if N0_on == True:
                if PrintStats==True:
                    print("Calculating redshift factors and intensities for n=0,  # of points", supergrid0.shape[0])
                i_g0 = gfactorf(a,i,supergrid0,mask0,sign0,rs0,Sub_Kep,BetaR,BetaPhi)
                I_0 = np.nan*np.zeros(rs0.shape[0])
                I_0[mask0] = intensity(rs0[mask0],i_g0[mask0],Alpha_r,R_int0) 

            if N1_on == True:
                if PrintStats==True:
                    print("Calculating redshift factors and intensities for n=1,  # of points", supergrid1.shape[0])
                i_g1 = gfactorf(a,i,supergrid1,mask1,sign1,rs1,Sub_Kep,BetaR,BetaPhi)
                I_1 = np.nan*np.zeros(rs1.shape[0])
                I_1[mask1] = intensity(rs1[mask1],i_g1[mask1],Alpha_r,R_int0) 

            if N2_on == True:
                if PrintStats==True:
                    print("Calculating redshift factors and intensities for n=2,  # of points", supergrid2.shape[0])
                i_g2 = gfactorf(a,i,supergrid2,mask2,sign2,rs2,Sub_Kep,BetaR,BetaPhi)
                I_2 = np.nan*np.zeros(rs2.shape[0])
                I_2[mask2] = intensity(rs2[mask2],i_g2[mask2],Alpha_r,R_int0) 

            gf_fldr="SpectraData%s%s%s_a_%s_i_%s/"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i)
            os.makedirs(path+gf_fldr, exist_ok=True)
            
            filename=path+gf_fldr +"gfactors%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s.h5"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi)

            h5f = h5py.File(filename, 'w')

            if N0_on == True:
                h5f.create_dataset('gf0', data=i_g0)

            if N1_on == True:
                h5f.create_dataset('gf1', data=i_g1)

            if N1_on == True:
                h5f.create_dataset('gf2', data=i_g2)

            h5f.close()
            
            
            if PrintStats==True:
                print("File ",filename," created.")
            
            if R_int0 == rms(a):
                Rint02='isco'
                
            elif R_int0 == rh(a):
                Rint02='horizon'
                
            elif R_int0 == rph(a):
                Rint02='pco'    
                
            else:
                Rint02=R_int0
            
            filename=path+gf_fldr +"Intensities%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s_sig_%s_r0_%s.h5" %(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi,Alpha_r,Rint02)

            h5f = h5py.File(filename, 'w')

            if N0_on == True:
                h5f.create_dataset('Io0', data=I_0)

            if N1_on == True:
                h5f.create_dataset('Io1', data=I_1)

            if N1_on == True:
                h5f.create_dataset('Io2', data=I_2)

            h5f.close()
            
            if PrintStats==True:
                print("File ",filename," created.") 
                