function [k1s,k2s,new_rs,radial_angles,area_factor] = GKTH_find_radial_ks(p,layers,opts)

%%% Description
% Subsamples k-space radially based on proximity to the Fermi surface
% Takes radial lines from kx=ky=0 to the edge of the first BZ
% Samples more frequently near the Fermi surface
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
%   opts        :   Optional keyword arguements
%       base_space     : The minimum spacing between points along the
%                        radial direction. Value is a ratio of the number
%                        of points determined by proximity to Fermi
%                        surface. A value of 1 means 50% of the points will
%                        be given by the even spacing, 50% by the proximity
%                        to the Fermi surface
%       width          : How aggressively to sample near the Fermi surface.
%                        A smaller number corresponds to sharper sampling
%                        (may miss features away from Fermi surface) and a
%                        low number means broader sampling (may miss sharp
%                        features at Fermi surface)
%       just_use_layer : Only use a particular layer or layers for the
%                        subsampling. Useful if e.g. you only care about
%                        the anomalous Green function and the sampling is
%                        skewed by a ferromagnet layer
%
%%% Outputs:
%   k1s         :   Array of kx points
%   k2s         :   Array of ky points
%   new_rs      :   Array of radial distances
%  radial_angles:   List of angles of the radial lines
%   area_factor :   Array of areas of k-space represented by each k-point

arguments
    p
    layers
    opts.base_space=0.1;
    opts.width=0.1;
    opts.just_use_layer=0;
end

if opts.just_use_layer~=0
    nlayers=1;
    layers=layers(opts.just_use_layer);
else
    nlayers=length(layers);
end

temp=num2cell(zeros(nlayers,1));
[layers.Delta_0]=temp{:};

rs=zeros(p.ntest,p.nradials);
eigenvalues=zeros(4*nlayers,p.ntest,p.nradials);

if p.lattice_symmetry == '4mm'
    radial_angles=linspace(0,(pi/4)+1e-5,p.nradials);
elseif p.lattice_symmetry == 'mm'
    radial_angles=linspace(0,(pi/2)+1e-5,p.nradials);
elseif p.lattice_symmetry == 'm'
    radial_angles=linspace(p.m_symmetry_line,p.m_symmetry_line+pi,p.nradials);
else
    radial_angles=linspace(0,2*pi,p.nradials);
end

eff_angles=mod(radial_angles,pi/2);
eff_angles(eff_angles>pi/4)=eff_angles(eff_angles>pi/4)-pi/2;
    
for i=1:p.nradials
    rs(:,i)=linspace(0,0.5/cos(eff_angles(i)),p.ntest);
end


k1s=2*pi/(p.a)*rs.*cos(radial_angles);
k2s=2*pi/(p.a)*rs.*sin(radial_angles);
base_m=GKTH_hamiltonian_k(p,k1s,k2s,layers);
for i=1:p.ntest
    for j=1:p.nradials
        eigenvalues(:,i,j)=eig(base_m(:,:,i,j));
    end
end

avg_spectrum = squeeze(prod(abs(eigenvalues),1).^(1/(4*nlayers)));

weights = exp(-abs(avg_spectrum)./opts.width) .* rs .* (gradient(rs')');
constant = sum(weights)/p.ntest * opts.base_space;
weights = constant + weights;

cumdist=cumsum(weights);

new_rs=zeros(p.nfinal,p.nradials);
for j=1:p.nradials
    new_rs(:,j)=interp1(cumdist(:,j),rs(:,j),linspace(min(cumdist(:,j)),max(cumdist(:,j)),p.nfinal));
end

area_factor = 2 * pi/(p.nradials-1) .* new_rs .* (gradient(new_rs')');% * (pi/(p.a))^2;
area_factor(:,1)=area_factor(:,1)/2;
area_factor(:,end)=area_factor(:,end)/2;
area_factor(1,:)=area_factor(1,:)/2;
area_factor(end,:)=area_factor(end,:)/2;

k1s=2*pi/(p.a)*new_rs.*cos(radial_angles);
k2s=2*pi/(p.a)*new_rs.*sin(radial_angles);