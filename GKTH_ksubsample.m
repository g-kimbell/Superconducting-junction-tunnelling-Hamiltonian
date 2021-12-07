function [density_grid,compute_grid] = GKTH_ksubsample(p,layers)

%%% Description
% Creates matrices used for grid-based k-space subsampling.
% Currently it does this by assuming the contribution to the Greens
% function roughly follows 1/xi, i.e. the k1,k2 which are close to the
% Fermi level for the different layers in the structure. If there are lots
% of layers with different Fermi surfaces this will not be very effective.
%
%%% Inputs:
%   p               :   Global Parameter object defining the stack and 
%                       calculation paramters.
%   layers          :   An array of Layer objects, defining the junction 
%                       structure.
%   point_fraction  :   The fraction of the number of points to aim for
%                       when subsampling. 1 is no subsampling, 0.1 is 10%
%                       the number of points. Default 0.1.
%   spread          :   A measure of how spread out or focussed the
%                       subsampling should be. Lower values e.g. 0.8 lead 
%                       to more aggressive high sampling/low sampling
%                       areas. High values e.g. 1.2 leads to larger regions
%                       of intermediate sampling rate. Default = 1.
%
%%% Outputs:
%   density_grid    :   An nxn matrix (defined by p.nkpoints), the value m
%                       in each element describes that it is part of an mxm
%                       block to calculate. So if m=1, the Greens function
%                       is calculated at every point. If m=2, the Greens
%                       function is calculate once per 2x2 set of points.
%   compute_grid    :   An nxn matrix (defined by p.nkpoints) of booleans.
%                       The value is 1 if it is in the bottom left corner
%                       of a sub-block defined by density_grid and 
%                       indicates that the Greens function should be 
%                       calculated here. A value of 0 indicates that the 
%                       Greens function calculation should be skipped.

spread = 1;
% 'spread' determines how spread out or focussed the subsampling should be.
% Lower values e.g. 0.8 lead to more aggressive high sampling/low sampling
% areas. High values e.g. 1.2 leads to larger regions of intermediate 
% sampling rate. I think it's probably best to leave it at 1.
block_sizes=[1,2,3,4]; %Final blcoks will be 2^block_sizes
%%% Function for making and refining the compute and density grids
function [npoints,density_grid,compute_grid] = make_grid(initial_spectrum,scale,spread)
    density_grid_inital=round(log((scale.*initial_spectrum).^(1/spread)));
    density_grid_inital(density_grid_inital<0)=0;
    density_grid_inital(density_grid_inital>max(block_sizes))=max(block_sizes);
    density_grid_initial=2.^density_grid_inital;
    % Refining the grid so it is in blocks of 1x1, 2x2, 4x4, 8x8, 16x16
    density_grid = ones([p.nkpoints,p.nkpoints]);
    compute_grid = ones([p.nkpoints,p.nkpoints]);
    for b = 2.^block_sizes
        for i = 1:p.nkpoints/b
            for j = 1:p.nkpoints/b
                if all(density_grid_initial(i*b-(b-1):i*b,j*b-(b-1):j*b)>=b)
                    density_grid(i*b-(b-1):i*b,j*b-(b-1):j*b)=b;
                    compute_grid(i*b-(b-1):i*b,j*b-(b-1):j*b)=0;
                    compute_grid(i*b-(b-1),j*b-(b-1))=1;
                end
            end
        end
    end
    if p.lattice_symmetry=='4mm'
        for i = 2:p.nkpoints
            for j = 1:i-1
                compute_grid(i,j)=0;
            end
        end
    elseif p.lattice_symmetry=='mm'
        % Do nothing
    elseif p.lattice_symmetry=='m'
        compute_grid(1:p.nkpoints/2,1:p.nkpoints/2)=0;
        compute_grid(p.nkpoints/2+1:end,p.nkpoints/2+1:end)=0;
    end
    npoints=sum(density_grid.^-2,'all');
end

nlayers=length(layers);
[eigenvalues] = GKTH_find_spectrum(p,layers);
if p.use_kspace_subsampling == true
    avg_spectrum =  prod(abs(eigenvalues),3).^(1/(4*nlayers));
    % Change scale until desired point fraction is achieved
    anonymous_function=@(x) make_grid(avg_spectrum,x,spread)/(p.nkpoints^2) - p.subsampling_point_fraction;
    options = optimset('TolX',0.01,'Tolfun',0.01); %'Display','iter',
    try
        best_scale = fzero(anonymous_function,[0,100],options);
    catch
        try
            best_scale = fzero(anonymous_function,[-100,1000],options);
        catch
            best_scale = fzero(anonymous_function,[0],options);
        end
    end
    [~,density_grid,compute_grid] = make_grid(avg_spectrum,best_scale,spread);
else
    density_grid = ones(p.nkpoints,p.nkpoints);
    compute_grid = ones([p.nkpoints,p.nkpoints]);
end

    
end


