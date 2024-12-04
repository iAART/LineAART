from aart_func import *
from params import * 

def weightf(xy,Is,lim, Nrho, Nvarphi):
    """
    Calculate the weight of each pixel in polar screen grid contributes to the line profle   
    :param xy: array of all {x,y} coordinates
    :param Is: array of all intesnity values
    :param lim: largest radius of on screen
    :param Nrho: number of rho screen coordinates
    :param Nvarphi: number of phi screen coordinates
    """
    rho = np.sqrt(xy[:,0]**2 + xy[:,1]**2)
    rhoC=rho_comp(rho)
    
    drhoC=rho_comp(lim)/(Nrho)
    drho=rho_by_rhoC(rhoC)*drhoC
    
    drho_dvarphi = drho * (2 * np.pi / Nvarphi)
    weight = Is * rho * drho_dvarphi
    
    return weight

def lp_hist(gs, weight, E_ref, binwidth):
    """
    Creates line profile histogram given pixel weights
    :param weight: weight of pixel
    :param E_ref: emission energy in rest frame
    :param binwidth: width of energy bins 
    """ 
    Emin=0
    Emax= E_ref*2
    energies=E_ref*gs
    
    binsize = np.arange(Emin, Emax, binwidth)
    
    # array of bin heights and bin boundaries
    heights, bin_boundaries = np.histogram(energies, bins=binsize, weights=weight)
    
    # midpoint of each bin
    bins = (bin_boundaries[1:] + bin_boundaries[:-1]) / 2
    
    #rescale heights so line profiles are rohbust against energy bin width
    heights=heights/binwidth
    
    return bins, heights

def line_profile(xy,gs,I,lim, Nrho, Nvarphi, E_ref, binwidth):
    weight = weightf(xy,I,lim, Nrho, Nvarphi)
    bins, heights = lp_hist(gs,weight, E_ref, binwidth)
    return bins, heights

def LineProfile(a,i,Sub_Kep,BetaR,BetaPhi,SigmaR,Rint0,Eref,rmin,rmax,binwidth):
    """
    Calculate the redshift factors and intensities factor from LensingBands, Rays, and gFactors files
    :param a: black hole spin
    :param i: observer inclination
    :params Sub_Kep, BetaR, BetaPhi: disk motion parameters

    :return: redshift factor and intensity at each point on on bserver screen.

    """
    if PrintStats==True:
        print("Computing Line Profile")


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
                lim0=int(h5f["lim0"][0])
                Nrho0=int(h5f["Nrho0"][0]) 
                Nvarphi0=int(h5f["Nvarphi0"][0]) 

            if N1_on == True:
                supergrid1=h5f['grid1'][:]
                mask1=h5f['mask1'][:]
                lim1=int(h5f["lim1"][0])
                Nrho1=int(h5f["Nrho1"][0])
                Nvarphi1=int(h5f["Nvarphi2"][0])               

            if N2_on == True:    
                supergrid2=h5f['grid2'][:]
                mask2=h5f['mask2'][:]
                lim2=int(h5f["lim2"][0])
                Nrho2=int(h5f["Nrho2"][0])
                Nvarphi2=int(h5f["Nvarphi2"][0])

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
            
            gf_fldr="SpectraData%s%s%s_a_%s_i_%s/"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i)
            
            fn=path+gf_fldr+"gfactors%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s.h5"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi)

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
           
            fn=path+gf_fldr+"Intensities%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s_sig_%s_r0_%s.h5"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi,SigmaR,Rint0)

            if PrintStats==True:
                print("Reading file: ",fn)

            h5f = h5py.File(fn,'r')

            if N0_on == True:
                I0=h5f['Io0'][:]

            if N1_on == True:
                I1=h5f['Io1'][:]

            if N2_on == True:    
                I2=h5f['Io2'][:]

            h5f.close() 
            
            if rmin == rms(a):
                rmin2='isco'
                
            elif rmin == rh(a):
                rmin2='horizon'
                
            elif rmin == rph(a):
                rmin2='pco'
                
            else:
                rmin2=rmin
            
            
            if rmax == rms(a):
                rmax2='isco'
                
            elif rmax == rh(a):
                rmax2='horizon'
                
            elif rmax == rph(a):
                rmax2='pco'    
                
            else:
                rmax2=rmax
                
                
            if Rint0 == rms(a):
                Rint02='isco'
                
            elif Rint0 == rh(a):
                Rint02='horizon'
                
            elif Rint0 == rph(a):
                Rint02='pco'    
                
            else:
                Rint02=Rint0
            
            filename=path+gf_fldr+"LineProfile%s%s%s_a_%s_i_%s_subkep_%s_Br_%s_Bphi_%s_sig_%s_r0_%s_Es_%s_rmin_%s_rmax_%s.h5"%(N_on_g[0],N_on_g[1],N_on_g[2],a,i,Sub_Kep,BetaR,BetaPhi,SigmaR,Rint02,Eref,rmin2,rmax2)
            h5f = h5py.File(filename, 'w')
            
            
            if N0_on == True:
                cond=(mask0) & (rs0>=rmin) & (rs0<=rmax)
                i_lp0x,i_lp0y=line_profile(supergrid0[cond],gf0[cond],I0[cond],lim0, Nrho0, Nvarphi0, Eref, binwidth)
                h5f.create_dataset('lp0', data=[i_lp0x,i_lp0y])
                
            if N1_on == True:
                cond=(mask1) & (rs1>=rmin) & (rs1<=rmax)
                i_lp1x,i_lp1y=line_profile(supergrid1[cond],gf1[cond],I1[cond],lim1, Nrho1, Nvarphi1, Eref, binwidth)
                h5f.create_dataset('lp1', data=[i_lp1x,i_lp1y])

            if N2_on == True: 
                cond=(mask2) & (rs2>=rmin) & (rs2<=rmax)
                i_lp2x,i_lp2y=line_profile(supergrid2[cond],gf2[cond],I2[cond],lim2, Nrho2, Nvarphi2, Eref, binwidth)
                h5f.create_dataset('lp2', data=[i_lp2x,i_lp2y])          
      
            h5f.create_dataset('binwidth', data=np.array([binwidth]))
            h5f.close()
            if PrintStats==True:
                print("File ",filename," created.") 
