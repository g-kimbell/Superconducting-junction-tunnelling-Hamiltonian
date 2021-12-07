function [lambda,layer] = GKTH_fix_lambda(p,layer,Delta_target)

%%% Description
% Finds and sets the BCS coupling constant to give a particular gap size
% at low temperature (0.1 K by default)
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
%  Delta_target :   The size of the gap at 0.1 K
%
%%% Outputs:
%
%   lambda      :   The BCS coupling constant found
%   layer       :   The layer object with updated BCS coupling constant
%

    if isa(layer,'Layer')
        layer.Delta_0=Delta_target;
        p.T=0.1*kB;
        p.h=0;
        [Fs_sums,~] = GKTH_Greens_radial(p,layer);
        lambda = Delta_target/(p.T*abs(Fs_sums));
        [layer.lambda]=lambda;
    else
        error("layer must be a single Layer object")
    end
end