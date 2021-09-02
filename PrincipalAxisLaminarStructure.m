function [n_parallel,n_perp]=PrincipalAxisLaminarStructure(n1,n2,pitch,DC, wavelength)
% [n_parallel,n_perp]=PrincipalAxisLaminarStructure(n1,n2,pitch,DC,wavelength)
%
% PrincipalAxisLaminarStructure solves plane wave propagation in an infinite
% periodic laminar structure.
% 
% This code can be found in the supplementary material of the review article:
% "A REVIEW OF SILICON SUBWAVELENGTH GRATINGS: BUILDING BREAK-THROUGH DEVICES WITH ANISOTROPIC METAMATERIALS", Nanophotonics.
%
% or 
%
% https://github.com/Photonics-RF-UMA/laminar_periodic_structure_homogenization
%
% Input parameters:
%   wavelength: free-space wavelength.
%   pitch     : period length in the same units as the wavelength.
%   DC        : duty cycle, which indicates the proportion of n1 in a period.
%   n1        : refractive index of the slab that has a thickness of DC*pitch. 
%   n2        : refractive index of the slab that has a thickness of (1-DC)*pitch. 
% For a graphical description of the parameters, see the figure below:
% 
% |          |             |          |             |                   
% |   n1     |    n2       |   n1     |    n2       |                   
% |          |             |          |             |         x                         
% |          |             |          |             |         ^          
% |          |             |          |             |         |          
% |          |             |          |             |         |          
% |          |             |          |             |  ...  y o----> z          
% |          |             |          |             |                   
% |          |             |          |             |                                                   
% |          |             |          |             |                   
% |          |             |          |             |                   
% | DC*pitch |(1-DC)*pitch |          |             |                   
% |<-------->|<----------->|          |             |                   
% |          |             |          |             |   
% 
% This function analytically solves the propagation of plane waves in a
% periodic laminar structure. The laminar structure is composed of periodically
% stacking slabs with refractive indices n1 and n2 along the z direction. The
% period is defined by the parameter *pitch* and the thicknesses of the
% slabs, which are DC*pitch and (1-DC)*pitch for the slabs with refractive 
% indices n1 and n2, respectively. See the schematic above for the geometry.
% 
% Specifically, this function computes the effective refractive indices that are seen by a plane wave 
% propagating along the principal axis of the infinite periodic laminar structure.
% The effective refractive indices that are computed are:
% 
%   n_parallel: the effective refractive index seen by a plane wave polarized
%               in x/y direction and propagating along the z direction.  
%   n_perp    : the effective refractive index seen by a plane wave polarized
%               in z direction and propagating along the x direction.
% 
% These two indices define an uniaxial anisotropic material in which the E
% and D field are related as follows:
% 
%        D=(Dx,Dy,Dz)=(n_parallel.^2*Ex+n_parallel.^2*Ey+n_perp.^2*Ez)*eps_0      
% 
% 
%  
%  (c) 2021 The authors
% 

a=pitch*DC;
b=pitch*(1-DC);

k0=2*pi/wavelength;

k1=n1*k0;
k2=n2*k0;

k1_z=@(Kx) sqrt(k1.^2-Kx.^2);
k2_z=@(Kx) sqrt(k2.^2-Kx.^2);

% delta parameter for TE polarization
delta_1=@(Kx) (n2.^2./n1.^2).*(k1_z(Kx)./k2_z(Kx));
delta_2=@(Kx) (n1.^2./n2.^2).*(k2_z(Kx)./k1_z(Kx));
delta=@(Kx) (delta_1(Kx)+delta_2(Kx))/2;

cos_kz_pitch=@(Kx) cos(k1_z(Kx).*a).*cos(k2_z(Kx).*b)-delta(Kx).*sin(k1_z(Kx).*a).*sin(k2_z(Kx).*b);

Kz=@(Kx) acos(cos_kz_pitch(Kx))./pitch;

% Kz when Kx=0
Kz_max=Kz(0);

%% Kx when Kz ==0

% We use as guess the Rytov formula 
Kx_max_rytov=(DC.*n1^-2+(1-DC)*n2.^-2)^-(1/2)*k0;

options = optimoptions('fsolve','Display','none');

% we solve for  cos_kz_pitch(Kx)=1 because this equation is better
% conditioned than the original Kx=0(ie acos(cos_kz_pitch(Kx))=0).
Kx_max=fsolve(@(Kx) 1-cos_kz_pitch(Kx),Kx_max_rytov,options);

%% final results
n_parallel= Kz_max/k0;
n_perp    = Kx_max/k0;

end

