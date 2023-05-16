Overview:

In this project, the transport of light into a silicon tube filled
with blood is simulated using Monte Carlo simulation. First, different 
geometries such as the wall thickness of the tube, the distance between 
the LED and the photodetector, and the wavelength of the photon are investigated.
In a next step, simulations are performed for different blood parameters 
such as oxygen saturation, blood volume fraction and hematocrit. 
The results of all simulations, which summarize the travel paths of the photons in the tube, 
the distribution of the luminous flux in the tube and the number of detected
photons, etc., are to be found in the results folder.

Brief instructions:

1. In MATLAB, run maketissue.m, which will create skinvessel_T.bin.
Any change of the geometrie parameter should be done here.
For chan
2. Compile mcxyz.c
3. Run mcxyz.c: mcxyz skinvessel . This will load the skinvessel_T.bin and skinvessel_H.mci 
needed by mcxyz.c, and run it.
4. In MATLAB, run lookmcxyz.m, which will read skinvessel_T.bin and skinvessel_F.bin and draw figures.

call tissue = makeTissueList(nm) : where nm specifies the wavelength of interest in [nm], 
will yield a data structure called tissue with the desired structure.

spectralLIB.mat : This file contains the optical absorption properties for blood 
and silicon, μa [cm-1], as functions of wavelength from 300 to 1000 nm.

mcxyz.c: is a computer simulation of the distribution of light 
within a complex tissue that consists of many different types of tissues, 
each with its own optical properties. The program uses the Monte Carlo 
method of sampling probabilities for the stepsize of photon movement
between scattering events and for the angles (θ,ψ) of photon scattering.
All boundary conditions are matched. Does not yet handle mismatched
boundaries.

mcxyz.c is written in ANSI standard C, and is easily compiled by any 
standard compiler. It reads input files (myname_T.bin, myname_H.mci) 
created by a MATLAB program maketissue.m, using a standard format. Therefore,
other programs can also prepare the input files.

The 3D Monte Carlo generates an output file of relative fluence rate, 
F(y,x,z) [W/cm2 per W delivered] or [1/cm2]. The spatial distribution of
absorption is obtained by the product of the fluence rate and the absorption 
coefficient: A(y,x,z) [1/cm3] = F(y,x,z) x muav(T(y,x,z), where muav(i) [cm-1]
is the absorption coefficient of the ith tissue type (the v is for voxel).

