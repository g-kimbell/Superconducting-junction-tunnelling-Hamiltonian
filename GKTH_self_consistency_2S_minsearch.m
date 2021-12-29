function[Delta,layers,history] = GKTH_self_consistency_2S_minsearch(p,layers,layers_to_check)

%%% Description
% GKTH_self_consistency_2S uses self-consistency to find the Delta values 
% for two different superconductors in a stack. This uses a much rougher
% approach than for a single superconductor and is more prone to error, and
% will not catch first order transitions.
% 
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
% layers_to_check : An array of integers defining which layers in the
%                   stack should have Delta calculated.
%
%%% Outputs:
%   Delta       :   Double, gap calcalated after running self-consistency.
%   layers      :   The array of Layer objects with updated Delta values
%   history     :   Matrix containing the checked Delta values along with
%                   the change in Delta from the self-consistency
%                   calculation. Used for analysing algorithm progress.
%

arguments
    p
    layers
    layers_to_check = [1 2]
end
tol = p.abs_tolerance_self_consistency_1S;

% This is the function used in the find-zero. The residual is the change in
% Delta after one iteration of self-consistency.
function[residual] = GKTH_self_consistency_2S_residual(x)
    for i=layers_to_check
        layers(i).Delta_0=x(i);
    end
    [Fs_sums,~] = GKTH_Greens_radial(p,layers);
    Deltas_iterate=[layers.lambda].*p.T.*abs(Fs_sums)';
    abs_change=Deltas_iterate-x;
    residual=norm(abs_change)^2;
    history=[history;[x,residual]];
end

%%%%%%%%%%%%%%%%%%% Root finding %%%%%%%%%%%%%%%%
history = [];
options=optimset('TolX',tol,'TolFun',tol);
initial_guess=[layers(layers_to_check(1)).Delta_0,layers(layers_to_check(2)).Delta_0];
Delta = fminsearch(@GKTH_self_consistency_2S_residual,initial_guess,options);

% Set the Delta for the layers being checked
for j = layers_to_check
    layers(j).Delta_0=real(Delta(j));
end

disp("Soltuion Deltas("+p.T+", "+p.ts(1)+") = "+Delta+" eV")
end