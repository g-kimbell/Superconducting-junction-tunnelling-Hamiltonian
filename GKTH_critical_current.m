function [jc,phase] = GKTH_critical_current(p,layers,opts)
%%% Description
% Finds the maximum current vs phase for a given structure
%
%%% Inputs
%   p      :    Global Parameter object defining the stack and 
%               calculation paramters.
%   layers :    An array of Layer objects defining the stack
%   opts   :    Optional keyword arguements
%       layer_to_vary : which layer to change the phase, default last
%       initial_guess : initial estimate of critical phase
%       maxCalcs      : maximum calculations in current calculation
%
%%% Outputs
%   jc     :    the critical current
%   phase  :    the critical phase

arguments
    p
    layers
    opts.layer_to_vary=NaN;
    opts.initial_guess=1.5;
    opts.maxCalcs=500;
end

% If no layer to vary set, make it the last layer
if isnan(opts.layer_to_vary)
    opts.layer_to_vary=length(layers);
end

% The function to minimize to get jc. Returns -jc so a minimizer
% can be used to find the maximum jc
function [js]=jc_function(xs)
    tic
    js=zeros(length(xs),1);
    for i=1:length(xs)
        x=xs(i);
        layers_temp=layers;
        layers_temp(opts.layer_to_vary).phi=x;
        [j_t] = GKTH_Greens_current_radial(p,layers_temp,maxCalcs=opts.maxCalcs);
        js(i)=-j_t(1);
    end
    toc
end

% The minimum search
options = optimset('TolX',1e-3,'Tolfun',1e10);
[phase,jc] = fminsearch(@jc_function,opts.initial_guess,options);
phase=mod(phase,2*pi);
% Force the result between -pi and pi
if phase>pi
    phase=phase-2*pi;
end
jc=-jc;
end