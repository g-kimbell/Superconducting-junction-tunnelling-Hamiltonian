function[xis] = GKTH_spectrum(p,L,mu,type)

%%% Description
% GKTH_spectrum gives the electronic disperson based on the grid-based
% k-sampling defined by the Global Params object p
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   L           :   The Layer object for which we are calculting the
%                   dispersion
%   mu          :   float, the chemical potential
%   type        :   string, either 'tb' or 'para' for tight-binding or
%                   parabolic dispersion respectively
%
%%% Outputs:
%
%   xis         :   Matrix containing the electronic dispersion values at
%                   each k-point in the grid defined by p
%

arguments
    p
    L
    mu
    type='tb'
end
if strcmp(type,'tb')
    xis = L.tNN*(cos(p.k1*p.a)+cos(p.k2*p.a)) + L.tNNN*cos(p.k1*p.a).*cos(p.k2*p.a) + mu;
elseif strcmp(type,'para')
    hbar=6.582119569e-16;
    me=9.10938356e-31;
    xis = hbar^2*((p.k1*p.a).^2+(p.k2*p.a).^2)/(2*me)+mu;
else
    error("Spectrum type must be 'tb' for tight binding or 'para' for parabolic.")
end