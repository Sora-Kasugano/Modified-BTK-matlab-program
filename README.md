# Modified-BTK-matlab-program
The matlab program is coded by github user Sora-Kasugano.
All rights reserved.

This program is intended for the experimental fit of modified BTK theory. Related work stemmed from the paper: . The following codes are based from this paper.
 
Tips of this program:
1: ''V'' stands for a group of DC bias voltage data. By default, the unit of V is mV.
2: ''aa'' stands for a row vector, which contains the following fit parameters.
a= aa(1); b= aa(2); Delta= aa(3); Gama=aa(4); Z= aa(5); P= aa(6); npanel= aa(7); T= aa(8);

     a ,b          = lower and upper limits of the integral.
     Delta        = Superconducting energy gap (mV).
     Gama       = Inelastic scattering factor with its unit mV.
     Z              = Barrier factor
     P              = Spin Polarization
     T              = Experimental Temperature (K)
     npanel      = number of panels to use in the integration, with Total number of nodes = 2*panel + 1
The output form are [S,DI1], where S represents the DC bias voltage, and DI1 the normalized differential conductance.    
