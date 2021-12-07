function[Tcs] = GKTH_Tc_finder(p,layers,layers_to_check)

%%% Description
% Finds the critical temperature of superconductors in a stack of layers
% !!! This seems to break when you tunnel between two superconductors   !!!
% !!! This also breaks when you have both a first and second order      !!!
% !!! transition since this will only catch the second order transition !!!
% !!! It is better to find the critical field in that case              !!!
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
% layers_to_check:  Which layers to search for the Tc of
%
%%% Outputs:
%   Tcs         :   The critical temperatures of the layers

% The value of Delta which is considered Tc
Delta_0 = 1e-10;
Tcs=zeros(length(layers_to_check),1);
i=1;

% The function to find-zero to get Tc
function[rel_Delta_change] = GKTH_Tc_finder_function(T)
    temp=num2cell(zeros(length(layers))+Delta_0); % This is silly but Matlab demands it
    [layers.Delta_0]=temp{:};
    p.T=T;
    [Fs_sums,~] = GKTH_Greens(p,layers);
    Fs_sum=Fs_sums(layer_to_check);
    rel_Delta_change=(layers(layer_to_check).lambda*p.T*abs(Fs_sum)-Delta_0)/Delta_0;
end

for layer_to_check = layers_to_check
    % Findzero the solution
    tolerance_T = 1e-5;
    tolerance_Delta = 1e-9;
    options = optimset('TolX',tolerance_T,'Tolfun',tolerance_Delta);

    try
        Tc=fzero(@GKTH_Tc_finder_function,[0.01,0.05],options);
    catch
        try
            Tc=fzero(@GKTH_Tc_finder_function,[0.001,0.01],options);
        catch
            try
                Tc=fzero(@GKTH_Tc_finder_function,[0.0001,0.001],options);
            catch
                warning("Can't find solution above 0.0001 eV")
                Tc=0;
            end
        end
    end
    disp("Layer "+layer_to_check+" Tc = "+Tc+" eV")
    Tcs(i)=Tc;
    i=i+1;
end

end