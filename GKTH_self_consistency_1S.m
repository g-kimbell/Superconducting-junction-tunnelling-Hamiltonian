function[Delta,layers,residual_history] = GKTH_self_consistency_1S(p,layers,layers_to_check,catch_first_order)

%%% Description
% GKTH_self_consistency_1S finds Delta for a single superconducting layer
% or several identical superconducting layers
% 
%%% Inputs:
%   
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
% layers_to_check : An array of integers defining which layers in the
%                   stack should have Delta calculated.
% catch_first_order: Boolean, whether or not to search for the first-order
%                    transition if the initial find zero fails.
%
%%% Outputs:
%   Delta       :   Double, gap calcalated after running self-consistency.
%   layers      :   The array of Layer objects with updated Delta values
% residual_history: Matrix containing the checked Delta values along with
%                   the change in Delta from the self-consistency
%                   calculation. Used for analysing algorithm progress.
%

arguments
    p
    layers
    layers_to_check = 1
    catch_first_order = true
end
tol = p.abs_tolerance_self_consistency_1S;
min_Delta=tol;
max_Delta=layers(layers_to_check(1)).Delta_0;

residual_history = [];

% First get the density grid. This shouldn't change between function calls
%[density_grid,compute_grid] = GKTH_ksubsample(p,layers);

% This is the function used in the find-zero. The residual is the change in
% Delta after one iteration of self-consistency.
function[residual] = GKTH_self_consistency_1S_residual(Delta_0_fits)
    residual=zeros(length(Delta_0_fits),1);
    for k = 1:length(Delta_0_fits)
        Delta_0_fit=Delta_0_fits(k);
        for i=layers_to_check
            layers(i).Delta_0=Delta_0_fit;
        end
        %[Fs_sums,~] = GKTH_Greens(p,layers,density_grid=density_grid,compute_grid=compute_grid);
        [Fs_sums,~] = GKTH_Greens_radial(p,layers,just_use_layer=layers_to_check);
        idx=layers_to_check(1);
        residual(k)=(layers(idx).lambda*p.T*abs(Fs_sums(idx))-Delta_0_fit);
        residual_history=[residual_history;[Delta_0_fit,residual(k)]];
    end
end

% Newton Raphson function used to catch the first order transition if the
% find zero fails.
function[x,fx,dfdx] = newton_raphson(x,fx,tol,maxItr)
    prev_x=x;
    prev_fx=fx;
    x=x+fx;
    dfdx=0;
    itr=0;
    min_neg_x=max_Delta;
    while itr<maxItr
        fx=GKTH_self_consistency_1S_residual(x);
        if fx>0
            opt = optimset('TolX',tol,'Tolfun',1e-99);
            try
                x=fzero(@GKTH_self_consistency_1S_residual,[x,min_neg_x],opt);
            catch
                try
                    x=fzero(@GKTH_self_consistency_1S_residual,[x,1.2*min_neg_x],opt);
                catch
                    warning("Failed when trying to swap to fzero. Can happen if there's too much noise. Try increasing the number of kpoints or tolerance in Greens calculation.")
                    x=0;
                end
            end
            break
        elseif fx<0 && x<min_neg_x
            min_neg_x=x;
        end
        dfdx=(fx-prev_fx)/(x-prev_x);
        if abs(fx)<tol
            if dfdx<0
                break
            else
                opt = optimset('TolX',tol,'Tolfun',1e-99);
                try
                    x=fzero(@GKTH_self_consistency_1S_residual,[x+2*fx/dfdx,min_neg_x],opt);
                catch
                    try
                        x=fzero(@GKTH_self_consistency_1S_residual,[x+1.5*fx/dfdx,min_neg_x],opt);
                    catch
                        try
                            x=fzero(@GKTH_self_consistency_1S_residual,[x+1.1*fx/dfdx,min_neg_x],opt);
                        catch
                            warning("Giving up trying to find positive part");
                            x=0;
                        end
                    end
                end
                break
            end
        end
        prev_x=x;
        prev_fx=fx;
        x=prev_x-prev_fx/dfdx;
        
        if dfdx>0
            x=0;
            break
        end
        if x<min_Delta
            x=0;
            break
        end
        itr=itr+1;
        if itr==maxItr
            x=NaN;
        end
    end
end

%%%%%%%%%%%%%%%%%%% Root finding %%%%%%%%%%%%%%%%

% x=linspace(0,5e-4,10);
% fx=GKTH_self_consistency_1S_residual(x);
% plot(x,fx);
% keyboard

% First try an interval search from min_Delta to max_Delta
x1=max(min_Delta,max_Delta/100);
x2=max_Delta;
f1=GKTH_self_consistency_1S_residual(x1);
f2=GKTH_self_consistency_1S_residual(x2);
options = optimset('TolX',tol,'Tolfun',1e-99);
if f1>0 && f2<0 % If +- then we can interval fzero
    Delta=gk_fzero(@GKTH_self_consistency_1S_residual,[x1,x2,f1,f2],options);
elseif f2>0 % If ++ or -+ then need to increase max_Delta then interval fzero
    if f1<0 && catch_first_order==false
        Delta=0;
    else
        itr_max_search=0;
        while f2>0
            x1=x2;
            f1=f2;
            x2=x2*2;
            f2=GKTH_self_consistency_1S_residual(x2);
            itr_max_search=itr_max_search+1;
            if itr_max_search==10
                error("Can't find a negative change in Delta in self-consistency eq. Try increasing initial Delta_0 value.")
            end
        end
        Delta=gk_fzero(@GKTH_self_consistency_1S_residual,[x1,x2,f1,f2],options);
    end
elseif f1<0 && f2<0 % If -- then Newton-raphson
    if catch_first_order==true
        [Delta,~,~] = newton_raphson(max_Delta,f2,tol,20);
    else
        Delta=0;
    end
end

% Set the Delta for the layers being checked
for j = layers_to_check
    layers(j).Delta_0=Delta;
end

disp("Soltuion Delta("+p.T+", "+p.h+") = "+Delta+" eV")

end