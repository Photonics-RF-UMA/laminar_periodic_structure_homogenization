# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 12:27:07 2020

@author: ah
"""
from scipy.optimize import fsolve
import numpy as np


# auxiliary function to solve complex functions of complex variables.
# i.e This function emulates Matlab's fsolve behaviour
def __fsolve_complex(fun,x0,*args,**kwargs):
    x0=np.array([np.real(x0),np.imag(x0)])
    def fun_vec(x):
        y=fun(x[0]+1j*x[1]);
        return np.array([np.real(y),np.imag(y)])
    x=fsolve(fun_vec,x0,*args,**kwargs);
    return [x[0]+1j*x[1],];


def PrincipalAxisLaminarStructure(n1,n2,pitch,DC, wavelength):
    """
    (n_parallel,n_perp) = PrincipalAxisLaminarStructure(n1,n2,pitch,DC,wavelength)
|

    PrincipalAxisLaminarStructure solves plane wave propagation in an  infinite
    periodic laminar structure.
    This code can be found in the supplementary material of the review article:
    "A REVIEW OF SILICON SUBWAVELENGTH GRATINGS: BUILDING BREAK-THROUGH DEVICES WITH ANISOTROPIC METAMATERIALS", Nanophotonics.
 
    or 

    https://github.com/Photonics-RF-UMA/laminar_periodic_structure_homogenization

    Input parameters:
    
    :param n1: float
        refractive index of the slab that has a thickness  of DC*pitch
    :param n2: float
        refractive index of the slab that has a thickness  of (1-DC)*pitch
    :param pitch: float
         period length in the same units as the wavelength
    :param DC: float
        duty cycle, indicates the proportion of n1 in a period.
    :param wavelength: float
        free space wavelength
    :return: tuple(n_parallel,n_perp)

    For a graphical description of the parameters see the figure below::

        |          |             |          |             |
        |   n1     |    n2       |   n1     |    n2       |
        |          |             |          |             |         x
        |          |             |          |             |         ^
        |          |             |          |             |         |
        |          |             |          |             |         |
        |          |             |          |             |  ...  y o----> z
        |          |             |          |             |
        |          |             |          |             |
        |          |             |          |             |
        |          |             |          |             |
        | DC*pitch |(1-DC)*pitch |          |             |
        |<-------->|<----------->|          |             |
        |          |             |          |             |

    This function analytically solves the propagation of plane waves in a
    periodic laminar structure. The laminar structure is composed of periodically
    stacking slabs with refractive indices n1 and n2 along the z direction. The
    period is defined by the parameter *pitch* and the thickness of the
    slabs, which are DC*pitch and (1-DC)*pitch for the slab with refractive 
    indices n1 and n2, respectively. See the schematic above for the geometry.

    Specifically, this function computes the effective refractive indices that are seen by plane wave
    propagating along the principal axis of the infinite periodic laminar structure.
    The effective refractive indices that are computed are:

      n_parallel: the effective refractive index seen by a plane wave polarized
                  in x/y direction and propagating along the z direction
      n_perp    : the effective refractive index seen by a plane wave polarized
                  in z direction and propagating along the x direction

    This two indices define an uniaxial anisotropic material in which the E
    and D field are related as follows:

           D=(Dx,Dy,Dz)=(n_parallel.^2*Ex+n_parallel.^2*Ey+n_perp.^2*Ez)*eps_0



     (c) 2021 The authors
    """
    a=pitch*DC;
    b=pitch*(1-DC);
    
    k0=2*np.pi/wavelength;
    
    k1=n1*k0;
    k2=n2*k0;
    
    k1_z = lambda Kx: np.sqrt(k1**2-Kx**2,dtype=complex);
    k2_z = lambda Kx: np.sqrt(k2**2-Kx**2,dtype=complex);
    
    
    # Delta_TE as defned in eq.XX of [1]
    delta_1 = lambda Kx: (n2**2/n1**2)*(k1_z(Kx)/k2_z(Kx));
    delta_2 = lambda Kx: (n1**2/n2**2)*(k2_z(Kx)/k1_z(Kx));
    delta   = lambda Kx: (delta_1(Kx)+delta_2(Kx))/2;
    
    
    # np.cos(Kz*pitch)=f(k_x)
    f = lambda Kx: np.cos(k1_z(Kx)*a)*np.cos(k2_z(Kx)*b)-delta(Kx)*np.sin(k1_z(Kx)*a)*np.sin(k2_z(Kx)*b);
    
    Kz=lambda Kx: np.arccos(f(Kx))/pitch;
    
    # Kz when Kx=0
    Kz_max=Kz(0);
    
    # Kx when Kz =0
    
    # We use as guess the Rytov formula 
    Kx_max_rytov=(DC*n1**-2+(1-DC)*n2**-2)**-(1/2)*k0;

    
    # we find the Kx value tha makes Kz=0 , that is Kx0 such that f(Kx0)=1
    
    g= lambda Kx: 1-f(Kx)
    
    Kx_max=__fsolve_complex(g,Kx_max_rytov);
    
    n_parallel= Kz_max/k0;
    n_perp    = Kx_max[0]/k0;
   
    return (n_parallel,n_perp)

if(__name__=="__main__"):
    n_parallel,n_perp=PrincipalAxisLaminarStructure(n1=3.476,n2=1.45,pitch=0.2,DC=0.5, wavelength=1.55)
    print("n_parallel={:5.4f}, n_perp={:5.4f}".format(n_parallel,n_perp))

    n_parallel, n_perp = PrincipalAxisLaminarStructure(n1=3.476, n2=1.45, pitch=0.6, DC=0.5, wavelength=1.55)
    print("n_parallel={:5.4f}, n_perp={:5.4f}".format(n_parallel, n_perp))
