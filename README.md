 # LineAART
This repository is an extended version of AART, a full detailed description of which is at https://github.com/iAART/aart.

This code computes a high-resolution relativistically broadened line emission from equatorial disks around a Kerr black hole.

## Table of Contents

- [Installation](#installation)
- [Key Features](#key-features)
- [Parameters](#parameters)
- [Example](#examples)

_______
### Installation

To use this simulation, clone the repository and install the required packages:

```bash
git clone https://github.com/iAART/LineAART.git
cd LineAART
pip install -r requirements.txt
```

Ensure you have Python 3.x and the following libraries installed:

    numpy
    matplotlib
    h5py
    imageio[v2]
    scipy
    scikit-image
_______
### Key Features

* **LensingBands:** sets up a grid of points on the observer screen of a spacially resolved image of the disk. The program can calculate the n = 0, 1, and 2 images. The grid is in polar coordates on the observer screen. 
* **RayTrace:** calculates the source radus for each observer screen pixel. 
* **Redshifts:** calculates the redshift factor g=E_o/E_s (ratio of photon energy in observer and source frames) for each observer screen pixel. 
* **Intensities:** calculates the intensity for each pixel on the observer screen. 
* **LineProfile:** calculates line profiles as a histogram over the observer screen pixels.

_______
### Parameters

The following parameters can be adjusted to modify the simulation:

|Parameter| Description                                                       |Type|
| ------- | ----------------------------------------------------------------- | ------------- |
| a | BH spin | ~ |
| i | observer inclination | ~ | 
| sub_kep | sub-Kepleniarity param | disk motion param |
| betr | radial velocity mixing param | disk motion param |
| betaphi | angular velocity param | disk motion param |
| sigma | disk emissivity fall off power | disk emmision param |
| r_0 | disk emissivity cut off | disk emission param |
| E_s | energy of emission | disk emission param |
| r_min | minimum disk radius | disk geometry param |
| r_max | maximum disk radius | disk geometry param |
| E_binwidth | width of line profiles bins | line profile param |
| limits | maximum radial screen coordinate | observer screen param |
| drho0, drho1, drho2 | compactified radial screen coordinate resolution for each lensing band n | observer screen param |
| dvarphi0, dvarphi1, dvarphi2 | angular screen coordinates resolution for each lensing band n | observer screen param |
| on_l0, on_l1, on_l2 | image number n on which to calculate lensing bands, on/off = 0 or 1 | ~ |
| on_r0, on_r1, on_r2 | image number n on which to calculate source radii, on/off = 0 or 1 | ~ |
| on_l0, on_l1, on_l2 | image number n on which to calculate redshift factors, intensities, and line profiles (on/off = 0 or 1) | ~ |
| path | the location of the saved data | ~ |
| comments | LineAART print out comments, on/off = 0 or 1 | ~ |

Refer to the code comments for detailed descriptions of each parameter and its units.

_______
### Examples

And example run of the code can be found as "LineAART.ipynb"
