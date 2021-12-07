function[Delta,layers,history] = GKTH_self_consistency_2S(p,layers,layers_to_check)

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

function[residual2D] = GKTH_self_consistency_2S_residual2D(x)
    for i=layers_to_check
        layers(i).Delta_0=x(i);
    end
    [Fs_sums,~] = GKTH_Greens_radial(p,layers);
    Deltas_iterate=[layers.lambda].*p.T.*abs(Fs_sums)';
    residual2D=real(Deltas_iterate-x);
end

%%%%%%%%%%%%%%%%%%% Root finding %%%%%%%%%%%%%%%%
% history = [];
% % First try an interval search from min_Delta to max_Delta
% options=optimset('TolX',tol,'TolFun',tol);
% initial_guess=[layers(layers_to_check(1)).Delta_0,layers(layers_to_check(2)).Delta_0];
% Delta = fminsearch(@GKTHS_self_consistency_2S_residual,initial_guess,options);

%%%%%%%%%%%%%%%%%% Newton Raphson %%%%%%%%%%%%%%%
history=[];
maxItr=100;
itr=0;
for i=1:length(layers_to_check)
    prev_x(i)=layers(layers_to_check(i)).Delta_0;
end
%prev_x=[layers(layers_to_check(1)).Delta_0,layers(layers_to_check(2)).Delta_0];
prev_fx=GKTH_self_consistency_2S_residual2D(prev_x);
x=prev_x+prev_fx;
fx=[1000 1000];
% xs=[prev_x , prev_fx];
while norm([fx,prev_fx])>tol && itr<maxItr
    fx=GKTH_self_consistency_2S_residual2D(x);
%     xs=[xs;x,fx];
    dfdx=(fx-prev_fx)./(x-prev_x);
    prev_x=x;
    prev_fx=fx;
    x=prev_x-prev_fx./dfdx;
    itr=itr+1;
end
Delta=x;


% Set the Delta for the layers being checked
for j = layers_to_check
    layers(j).Delta_0=real(Delta(j));
end

%disp("Soltuion Deltas("+p.T+", "+p.h+") = "+Delta+" eV")
end