import sys
import os
import argparse
import fileinput

parent_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
params_file = os.path.join(parent_directory, 'params.py')

def nullable_string(val):
    if not val:
        return "None"
    return val

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='AART')

    parser.add_argument('--a', help='BH Spin e.g., 0.94')

    parser.add_argument('--i', help="Observers' inclination [degrees] e.g., 17")



    parser.add_argument('--subkep', help="Sub-Keplerianity")

    parser.add_argument('--betar', help="Radial velocity factor")

    parser.add_argument('--betaphi', help="Angular velocity factor")  

    parser.add_argument('--sigma', help="Source emissitivity fall off")  

    parser.add_argument('--r_0', help="Source emissitivity zero")

    parser.add_argument('--r_min', help="min disk radius for spectrum")

    parser.add_argument('--r_max', help="max disk radius for spectrum")  

    parser.add_argument('--E_s', help="emission energy in disk's particle rest frame")

    parser.add_argument('--E_binwidth', help="energy bin width for spectrum")  

    
    
    parser.add_argument('--limits', help="Screen radius")  

    parser.add_argument('--drho0', help="Screen resolution, radial, n=0")

    parser.add_argument('--drho1', help="Screen resolution, radial, n=1")

    parser.add_argument('--drho2', help="Screen resolution, radial, n=2")

    parser.add_argument('--dvarphi0', help="Screen resolution, angular, n=0")

    parser.add_argument('--dvarphi1', help="Screen resolution angular, n=1")

    parser.add_argument('--dvarphi2', help="Screen resolution angular, n=2")
    
    
    
    parser.add_argument('--on_l0', help="n=0 lensingband on/off")

    parser.add_argument('--on_l1', help="n=1 lensingband on/off")

    parser.add_argument('--on_l2', help="n=2 lensingband on/off")

    parser.add_argument('--on_r0', help="n=0 source radius on/off")

    parser.add_argument('--on_r1', help="n=1 source radius on/off")

    parser.add_argument('--on_r2', help="n=2 source radius on/off")

    parser.add_argument('--on_g0', help="n=0 redshift factors on/off")

    parser.add_argument('--on_g1', help="n=1 redshift factors on/off")

    parser.add_argument('--on_g2', help="n=2 redshift factors on/off")



    parser.add_argument('--path',  help="Path to folder")

    parser.add_argument('--comments',  help="Print out comments")




########################################################

#   parser.add_argument('--n', help="Lensing band")

#   parser.add_argument('--maxbl', help="Max Baseline")

#   parser.add_argument('--gamma', help="Astro param gamma")
#   parser.add_argument('--mu', help="Astro param mu")
#   parser.add_argument('--sigma', help="Astro param sigma")

#   parser.add_argument('--armangle', help="Anisotropy direction")

#   parser.add_argument('--noise', help="Fluctuation scale")

#   parser.add_argument('--xyres', help="x and y number of points")

#   parser.add_argument('--trs', help="time number of points")

#   parser.add_argument('--xcorr', help="Spatial correlation in the x direction")

#   parser.add_argument('--xycorr', help="Spatial correlation in the y direction")

#   parser.add_argument('--snapshots', help="Number of Snapshots")

#   parser.add_argument('--production', type=int, default=1, help="Are we doing science?")

########################################################

    args = parser.parse_args()

    for line in fileinput.input(params_file, inplace=True):

        if line.strip().startswith('spin_case=') and args.a!=None:
            line = 'spin_case=%s\n'%args.a

        if line.strip().startswith('i_case=') and args.i!=None:
            line = 'i_case=%s\n'%args.i



        if line.strip().startswith('sub_kep=') and args.subkep!=None:
            line = 'sub_kep=%s\n'%args.subkep

        if line.strip().startswith('betar=') and args.betar!=None:
            line = 'betar=%s\n'%args.betar

        if line.strip().startswith('betaphi=') and args.betaphi!=None:
            line = 'betaphi=%s\n'%args.betaphi

        if line.strip().startswith('sigma=') and args.sigma!=None:
            line = 'sigma=%s\n'%args.sigma

        if line.strip().startswith('r_0=') and args.r_0!=None:
            line = 'r_0=%s\n'%args.r_0

        if line.strip().startswith('r_min=') and args.r_min!=None:
            line = 'r_min=%s\n'%args.r_min

        if line.strip().startswith('r_max=') and args.r_max!=None:
            line = 'r_max=%s\n'%args.r_max

        if line.strip().startswith('E_s=') and args.E_s!=None:
            line = 'E_s=%s\n'%args.E_s    

        if line.strip().startswith('E_binwidth=') and args.E_binwidth!=None:
            line = 'E_binwidth=%s\n'%args.E_binwidth



        if line.strip().startswith('limits=') and args.limits!=None:
            line = 'limits=%s\n'%args.limits

        if line.strip().startswith('drho0=') and args.drho0!=None:
            line = 'drho0=%s\n'%args.drho0

        if line.strip().startswith('drho1=') and args.drho1!=None:
            line = 'drho1=%s\n'%args.drho1

        if line.strip().startswith('drho2=') and args.drho2!=None:
            line = 'drho2=%s\n'%args.drho2

        if line.strip().startswith('dvarphi0=') and args.dvarphi0!=None:
            line = 'dvarphi0=%s\n'%args.dvarphi0

        if line.strip().startswith('dvarphi1=') and args.dvarphi1!=None:
            line = 'dvarphi1=%s\n'%args.dvarphi1

        if line.strip().startswith('dvarphi2=') and args.dvarphi2!=None:
            line = 'dvarphi2=%s\n'%args.dvarphi2



        if line.strip().startswith('on_l0=') and args.on_l0!=None:
            line = 'on_l0=%s\n'%args.on_l0

        if line.strip().startswith('on_l1=') and args.on_l1!=None:
            line = 'on_l1=%s\n'%args.on_l1

        if line.strip().startswith('on_l2=') and args.on_l2!=None:
            line = 'on_l2=%s\n'%args.on_l2

        if line.strip().startswith('on_r0=') and args.on_r0!=None:
            line = 'on_r0=%s\n'%args.on_r0

        if line.strip().startswith('on_r1=') and args.on_r1!=None:
            line = 'on_r1=%s\n'%args.on_r1

        if line.strip().startswith('on_r2=') and args.on_r2!=None:
            line = 'on_r2=%s\n'%args.on_r2   

        if line.strip().startswith('on_g0=') and args.on_g0!=None:
            line = 'on_g0=%s\n'%args.on_g0

        if line.strip().startswith('on_g1=') and args.on_g1!=None:
            line = 'on_g1=%s\n'%args.on_g1

        if line.strip().startswith('on_g2=') and args.on_g2!=None:
            line = 'on_g2=%s\n'%args.on_g2



        if line.strip().startswith('path=') and args.path!=None:
            line = 'path="'+args.path+'"\n'

        if line.strip().startswith('PrintStats=') and args.comments!=None:
            line = 'PrintStats=%s\n'%args.comments


       
        
        sys.stdout.write(line)
    
print("params.py updated")