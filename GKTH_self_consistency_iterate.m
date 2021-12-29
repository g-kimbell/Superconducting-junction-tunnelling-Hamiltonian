function[Deltas,layers,Deltas_itr,fs] = GKTH_self_consistency_iterate(p,layers)
%%% Description
% GKTH_self_consistency_iterate finds Delta of every layer in a stack by
% iteratively running the self-consistency equation.
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
%
%%% Outputs:
%   Deltas      :   Array containing the calculated gap of each layer in
%                   the stack
%   layers      :   The layers with modified Delta after self-consistency
%   Deltas_itr  :   Array containing the calculated gap of each layer in
%                   the stack at every iteration.
%   fs          :   1D Array of the resolution factor at each iteration
%

nlayers=length(layers);
maxItr=750;
iterate_factor=1; %Multiply the change by this each time. Might overshoot.
Deltas=[layers.Delta_0];
Deltas_itr=zeros(nlayers,maxItr);
n=0;
fs=zeros(1,maxItr);

for resolution_factor = [1]
    p_temp=p;
    p_temp.nfinal=p.nfinal/resolution_factor;
    p_temp.nradials=p.nradials/resolution_factor;
    rel_change = 1000;
    abs_change = 1000;
    rel_tol=p.rel_tolerance_self_consistency_iterate*10^(2*resolution_factor-2);
    abs_tol=p.abs_tolerance_self_consistency_iterate*10^(2*resolution_factor-2);
    while rel_change>rel_tol && norm(abs_change)>abs_tol && n<maxItr
        [Fs_sums,~] = GKTH_Greens_radial(p_temp,layers,layers_to_check=[1,2]);
        Deltas_iterate=[layers.lambda].*p_temp.T.*abs(Fs_sums)';
        abs_change=Deltas_iterate-Deltas;
        Deltas_new=Deltas+abs_change;
        rel_change=norm(Deltas_iterate./Deltas - 1);
        temp=num2cell(Deltas_new); % This is silly but Matlab demands it
        [layers.Delta_0]=temp{:};
        Deltas=Deltas_new;
        n=n+1;
        Deltas_itr(:,n)=Deltas; % In case you want to check the progress
        fs(n)=resolution_factor; % In case you want to check the progress
    end
end
if n==maxItr
    warning("Hit max itr")
end
end