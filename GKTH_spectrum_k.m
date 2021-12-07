function[xis] = GKTH_spectrum_k(p,L,k1s,k2s)

%%% Description
% GKTH_spectrum gives the electronic disperson for any k values defined by
% the k1s and k2s arrays
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   L           :   The Layer object for which we are calculting the
%                   dispersion. The layer object contains the chemical
%                   potential and type of dispersion
%
%%% Outputs:
%
%   xis         :   Matrix containing the electronic dispersion values at
%                   each k-point defined by k1s and k2s
%

arguments
    p
    L
    k1s
    k2s
end
if strcmp(L.dispersion_type,'tb')
    xis = L.tNN*(cos(k1s*p.a)+cos(k2s*p.a)) + L.tNNN*cos(k1s*p.a).*cos(k2s*p.a) + L.mu;
elseif strcmp(L.dispersion_type,'para')
    hbar=6.582119569e-16;
    me=9.10938356e-31;
    xis = hbar^2*((k1s*p.a).^2+(k2s*p.a).^2)/(2*me)+L.mu;
else
    error("Spectrum type must be 'tb' for tight binding or 'para' for parabolic.")
end