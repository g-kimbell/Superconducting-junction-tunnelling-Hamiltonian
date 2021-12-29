function[Delta,layers,history] = GKTH_self_consistency_2S_taketurns(p,layers,layers_to_check)

%%% Description
% GKTH_self_consistency_2S uses self-consistency to find the Delta values 
% for two different superconductors in a stack. This uses a much rougher
% approach than for a single superconductor and is more prone to error, and
% will not catch first order transitions.
% Here, we do a numerical Newton-Raphson to find the approximate solution
% to self-consistency of one superconductor whilst holding the other
% superconductor gap constant. We then repeat the process with the other
% superconductor and keep repeating until within a tolerance of the
% solution. This seems more robust than methods like gradient descent.
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

function[residual2D] = GKTH_self_consistency_2S_residual2D(x)
    for i=layers_to_check
        layers(i).Delta_0=x(i);
    end
    [Fs_sums,~] = GKTH_Greens_radial(p,layers,layers_to_check=[1,2]);
    Deltas_iterate=[layers.lambda].*p.T.*abs(Fs_sums)';
    residual2D=real(Deltas_iterate-x);
    history=[history;[x,residual2D]];
end

%%%%%%%%%%%%%%%%%% Newton Raphson %%%%%%%%%%%%%%%
history=[];
minItr=10;
maxItr=100;
itr=0;
for i=1:length(layers_to_check)
    x(i)=layers(layers_to_check(i)).Delta_0;
end

fx=[1000 1000];
xshift=[1000 1000];
dx=tol/2;
while (norm(xshift)>tol && itr<maxItr) || itr<minItr
    fx=GKTH_self_consistency_2S_residual2D(x);
	if mod(itr,2)
        dfx=GKTH_self_consistency_2S_residual2D(x+[dx,0]);
        dfdx=(dfx(1)-fx(1))/dx;
        x=x-[fx(1)/dfdx,0];
        xshift(1)=fx(1)/dfdx;
    else
        dfx=GKTH_self_consistency_2S_residual2D(x+[0,dx]);
        dfdx=(dfx(2)-fx(2))/dx;
        x=x-[0,fx(2)/dfdx];
        xshift(2)=fx(2)/dfdx;
    end
    itr=itr+1;
end
Delta=x;

% Set the Delta for the layers being checked
for j = layers_to_check
    layers(j).Delta_0=real(Delta(j));
end

end