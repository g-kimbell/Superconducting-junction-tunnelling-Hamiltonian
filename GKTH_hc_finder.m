function[hcFO,hcSO,values,all_history] = GKTH_hc_finder(p,layers,layers_to_check,max_Delta,catch_first_order)

%%% Description
%
%%% Inputs
%   p              :   Global Parameter object defining the stack and 
%                      calculation paramters
%   layers         :   An array of Layer objects, defining the junction 
%                      structure.
% layers_to_check  :   Which layers to adjust the Delta value of
%   max_Delta      :   A guess of the maximum possible Delta value
% catch_first_order:   bool, whether or not to try searching for a first 
%                      order transition
%
%%% Outputs
%   hcFO        :   The first-order critical field
%   hcSO        :   The second-order critical field
%   values      :   A 5 x n array containing the guessed Delta along with
%                   the maximum and minimum possible gap and field at each 
%                   step. Used to analyse progress of the function.
%   all_history :   The field and function values at every point checked.

arguments
    p
    layers
    layers_to_check
    max_Delta
    catch_first_order=true
end

%%%%%%%%%%%%%%%%% Initialise variables, define functions %%%%%%%%%%%%%%%%%%

% The value of Delta which is considered hc
Delta_0 = 1e-10;
p.abs_tolerance_Greens=1e-99; % need to set this low because Delta0 is so low
superconducting=false;

% The function to find-zero to get second-order hc
function[rel_Delta_change] = GKTH_SOhc_finder_function(h)
    for j=layers_to_check
        layers(j).Delta_0=Delta_0;
    end
    p.h=h;
    [Fs_sums] = GKTH_Greens_radial(p,layers);
    Fs_sum=Fs_sums(layers_to_check(1));
    rel_Delta_change=(layers(layers_to_check(1)).lambda*p.T*abs(Fs_sum)-Delta_0);
end

% The function to minimize to get first-order hc
function[rel_Delta_change] = GKTH_FOhc_finder_function(Delta_0s)
    rel_Delta_change=zeros(1,length(Delta_0s));
    for i=1:length(Delta_0s)
        for j=layers_to_check
            layers(j).Delta_0=Delta_0s(i);
        end
        [Fs_sums] = GKTH_Greens_radial(p,layers);
        Fs_sum=Fs_sums(layers_to_check(1));
        rel_Delta_change(i)=-(layers(layers_to_check(1)).lambda*p.T*abs(Fs_sum)-Delta_0s(i));
    end
end

% Function for terminating the search
function [stop] = outfun(x,optimvalues,~)
    history=[history; x,optimvalues.fval];
    if optimvalues.fval<0
        superconducting = true;
        stop = true;    
    else
        stop = false;
        if length(history)>3
            hisorty=sortrows(history);
            [min_dD,idx] = min(hisorty(:,2));
            if idx>1 && idx<length(history)
                dfdx=diff(hisorty(:,1))./diff(hisorty(:,2));
                if max(abs(dfdx(idx-1:idx)))*(hisorty(idx+1,2)-hisorty(idx-1,2))<min_dD
                    superconducting = false;
                    stop = true;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% SECOND ORDER TRANSITION %%%%%%%%%%%%%%%%%%%%%%%%%

tolerance_h = p.abs_tolerance_hc;
h_max = layers(layers_to_check(1)).Delta_0*2;
options = optimset('TolX',tolerance_h,'Tolfun',1e-99,'Display','iter');
try
    hcSO=fzero(@GKTH_SOhc_finder_function,[0,h_max],options);
catch     
    try
        hcSO=fzero(@GKTH_SOhc_finder_function,[0,10*h_max],options);
    catch
        hcSO=0;
    end
end
disp("SO hc(T="+p.T+") = "+hcSO+" eV")
hcFO=hcSO;

%%%%%%%%%%%%%%%%%%%%%%%%% FIRST ORDER TRANSITION %%%%%%%%%%%%%%%%%%%%%%%%%%

verbose=false;
values=[];
all_history=[];

tolerance_h=hcSO/1000; % Recalibrate the h tolerance

if hcSO ~= 0 && catch_first_order==true
    Deltas=[linspace((10*Delta_0)^0.5,(max_Delta)^0.5,10).^2, 2*max_Delta];
    dDs=zeros(1,length(Deltas));
    p.h=hcSO+tolerance_h; % to search very slightly above hSO
    [density_grid,compute_grid] = GKTH_ksubsample(p,layers);
    sign_changes=0;
    for k=1:length(Deltas)
        dDs(k)=GKTH_FOhc_finder_function(Deltas(k),density_grid,compute_grid);
        if sign_changes==0 && dDs(k)>0
            sign_changes=1;
        elseif sign_changes==1 && dDs(k)<0
            sign_changes=2;
        elseif sign_changes==2 && dDs(k)>0
            break
        end
    end
    signs=sign(dDs); % -1 is stable +1 is unstable
    if all(signs==1)
        % do nothing
    elseif signs(end)==-1
        error("Couldn't find a max Delta in FO transition search. Try increasing max Delta estimate")
    else
        max_Delta_idx = find(signs==-1,1,'last')+1;
        min_Delta_idx = find(signs(1:max_Delta_idx-1)==1,1,'last');
        max_Delta=Deltas(max_Delta_idx);
        if isempty(min_Delta_idx)
            min_Delta=0;
        else
            min_Delta=Deltas(min_Delta_idx);
        end
        
        if verbose == true
            hold on
            yline(hcFO)
            xline(min_Delta)
            xline(max_Delta)
            zs=zeros(k,1);
            all_history=[zs+hcSO+tolerance_h,Deltas(1:k)',dDs(1:k)',zs+hcSO,zs+inf,zs,zs+inf];
        end

        Delta_guess=(max_Delta+min_Delta)/2;
        Delta=Delta_guess;
        dx=1000;
        min_N_h=NaN;
        max_S_h=hcSO;
        h=1.05*hcSO;
        itr=0;
        if verbose==true
            values=[Delta_guess,min_Delta,max_Delta,max_S_h,min_N_h];
        end
        while dx>tolerance_h
            p.h=h;
            [density_grid,compute_grid] = GKTH_ksubsample(p,layers);
            options = optimset('MaxFunEvals',25,'display','off','OutputFcn', @outfun,'TolX',1e-8,'Tolfun',1e-8);
            guess_dD=GKTH_FOhc_finder_function(Delta_guess,density_grid,compute_grid);
            history=[Delta_guess,guess_dD];
            if guess_dD<0
                exitflag=-1;
                superconducting = true;
            else
                [Delta,~,exitflag] = fminbnd(@(x) GKTH_FOhc_finder_function(x,density_grid,compute_grid),min_Delta,max_Delta,options);
            end
            itr=itr+1;
            if verbose==true
                zs=zeros(size(history,1),1);
                all_history=[all_history;zs+h,history,zs+max_S_h,zs+min_N_h,zs+min_Delta,zs+max_Delta];
                if exitflag==-1 && superconducting == true
                    msg="superconducting.";
                else
                    msg="non-superconducting.";
                end
                disp("Iteration: "+itr)
                disp("Exit flag: "+exitflag)
                disp("Max superconducting h: "+max_S_h)
                disp("Min normal h: "+min_N_h)
                disp("h: "+h+" eV, this is "+msg)
                hold on
            end
            if exitflag==-1 && superconducting == true
                if isnan(min_N_h)
                    max_S_h=h;
                    h=1.05*h;
                else
                    max_S_h=h;
                    h=max_S_h+0.25*(min_N_h-max_S_h);
                end
                Delta_guess=Delta;
                if length(history)>1
                    N_Deltas=history(history(:,2)>0,1);
                    new_max_Delta=min(N_Deltas(N_Deltas > Delta));
                    new_min_Delta=max(N_Deltas(N_Deltas < Delta));
                    if ~isempty(new_max_Delta)
                        max_Delta=new_max_Delta;
                    end
                    if ~isempty(new_min_Delta)
                        min_Delta=new_min_Delta;
                    end
                end
            else
                min_N_h=h;
                h=max_S_h+0.25*(min_N_h-max_S_h);
            end
            hcFO = h;
            dx = min_N_h-max_S_h;
            if isnan(dx)
                dx=10000;
            end
            if verbose == true
                hold on
                try
                    yline(max_S_h)
                end
                try
                    yline(min_N_h)
                end
                xline(min_Delta)
                xline(max_Delta)
                values=[values ; Delta_guess,min_Delta,max_Delta,max_S_h,min_N_h];
            end
        end
    end
end
disp("FO hc(T="+p.T+") = "+hcFO+" eV")
end