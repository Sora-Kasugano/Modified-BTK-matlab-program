function  [S, DI1] = mBTK(aa,V)
% The matlab program is coded by github user Sora-Kasugano.
% All rights reserved.
%
% This program is intended for the experimental fit of modified BTK theory. Related work stemmed from the
% paper: . The following codes are based from this paper.
% 
% Tips of this program:
% 1: ''V'' stands for a group of DC bias voltage data. By default, the unit of V is mV.
% 2: ''aa'' stands for a row vector, which contains the following fit parameters.
% a= aa(1); b= aa(2); Delta= aa(3); Gama=aa(4); Z= aa(5); P= aa(6); npanel= aa(7); T= aa(8);
%
%      a ,b          = lower and upper limits of the integral.
%      Delta        = Superconducting energy gap (mV).
%      Gama       = Inelastic scattering factor with its unit mV.
%      Z              = Barrier factor
%      P              = Spin Polarization
%      T              = Experimental Temperature (K)
%      npanel      = number of panels to use in the integration, with
%                         Total number of nodes = 2*panel + 1
%
a= aa(1); b= aa(2); Delta= aa(3); Gama=aa(4); Z= aa(5); P= aa(6); npanel= aa(7); T= aa(8); 
% Every V gets a value of Andreev Reflection Current. 
for m=1:length(V)
    n = 2 * npanel + 1;                                          % total number of nodes
    h = (b-a)/(n-1);                                                % stepsize
    E = a:h:b;                                                         % divide the interval
    u02 = 0.5 * (1 + real(sqrt(((E - i * Gama).^2 - Delta^2)./(E - i * Gama).^2)));             % coefficient of BCS u0^2
    v02 = 0.5 * (1 - real(sqrt(((E - i * Gama).^2 - Delta^2)./(E - i * Gama).^2)));              % coefficient of BCS v0^2
    gama2 = (u02 + Z^2 * (u02 - v02)).^2;           % coefficient of BTK
    % define the coefficient of BTK theory A, B
    j = 1;
    while j <= length(E)
        if abs(E(j)-i*Gama) <= Delta
            AN(j) = real(Delta^2/((E(j) - i * Gama)^2 + (Delta^2 - (E(j) - i * Gama)^2) * (1 + 2 * Z^2)^2));
            BN(j) = 1 - AN(j);
            BP(j) = 1;
        else
            AN(j) = u02(j) * v02(j)/gama2(j);
            BN(j) = (u02(j) - v02(j))^2 * Z^2 * (1 + Z^2)/gama2(j);
            BP(j) = BN(j);
        end
        j=j + 1;
    end
    f = 1./(1+exp(11.5942*(E-V(m))/T)) - 1./(1+exp(11.5942*E/T));                                % f(E-eV)-f(E)
    z = ((1-P) * (1 + AN - BN) + P * (1 - BP)).* f;
    I(m) = (h/3)*(z(1)+4*sum(z(2:2:n-1))+2*sum(z(3:2:n-2))+z(n));
end
% Claculate the normalized differential conductance and output them with
% the form [S DI1], where S represents the DC bias voltage, and DI1 the
% normalized differential conductance.
S=V;
S(1)=[];
DI=diff(I);
DI1=DI/DI(length(V)-1);
% Since S and DI1 are in a same dimension, it's easy to plot the result fit of modified BTK theory.
