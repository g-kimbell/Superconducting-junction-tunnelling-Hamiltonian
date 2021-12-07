function[Ds] = GKTH_Delta_k(p,L,k1s,k2s)
%%% Description:
% GKTH_Delta builds the gap k-dependence, depending on the gap size and
% symmetry. Uses any grid of k-points supplied by k1s and k2s.
%
%%% Inputs:
%   p       :   Global Parameter object defining the stack and 
%               calculation paramters.
%   L       :   A Layer object, the layer the gap should be calculated for.
%   k1s     :   2D array of k points in x
%   k2s     :   2D array of k points in y
%
%%% Outputs:
%   Ds      :   An matrix of the gap with the same dimensions as k1s or k2s

if L.symmetry == "d"
    Ds = L.Delta_0/2 * (cos(k1s*p.a) - (cos(k2s*p.a)));
elseif L.symmetry == "s"
    Ds = zeros(size(k1s))+L.Delta_0;
elseif L.symmetry == "n"
    Ds = zeros(size(k1s));
else
    warning("Input a symmetry for the gap. s = s-wave, d = d-wave, n=non-superconducting. p-wave not supported yet.")
end