function[Ds] = GKTH_Delta(p,L,Delta)
%%% Description:
% GKTH_Delta builds the gap k-dependence, depending on the gap size and
% symmetry. Uses the square grid of k points defined in p.
%
%%% Inputs:
%   p       :   Global Parameter object defining the stack and 
%               calculation paramters.
%   L       :   A Layer object, the layer the gap should be calculated for.
%
%   Delta   :   Double, the size of the gap.
%
%%% Outputs:
%   Ds      :   An nxn matrix (dependent on nkpoints in p) with the size of
%               the gap at each k1, k2 pair

if L.symmetry == "d"
    Ds = Delta/2 * (cos(p.k1*p.a) - (cos(p.k2*p.a)));
elseif L.symmetry == "s"
    Ds = zeros(p.nkpoints,p.nkpoints)+Delta;
elseif L.symmetry == "n"
    Ds = zeros(p.nkpoints,p.nkpoints);
else
    warning("Input a symmetry for the gap. s = s-wave, d = d-wave, n=non-superconducting. p-wave not supported yet.")
end
