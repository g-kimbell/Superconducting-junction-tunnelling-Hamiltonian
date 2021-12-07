function [js,E_kresolved_final] = GKTH_Greens_current(p,layers,opts)
%%% Description
% GKTH_Greens_current inverts the Hamiltonian and takes a sum over k1, k2 
% and Matsubara frequencies, returning the current between each layer.
%
%%% Inputs:
%   p            :   Global Parameter object defining the stack and 
%                    calculation paramters
%   layers       :   An array of Layer objects, defining the junction 
%                    structure.
%
%%% Outputs:
%
%   jz_up        :   An array of the up-up current between each pair of
%                    layers in the stack.
%   jz_down      :   An array of the down-down current between each pair
%                    of layers in the stack.
%
%%% Function:
arguments
    p
    layers
    opts.include_spin = false;
    opts.density_grid = NaN;
    opts.compute_grid = NaN;
    opts.maxCalcs=500;
    opts.maxMatsubara=1e7; % before 1e6
    opts.verbose=false;
end
include_spin=opts.include_spin;
maxCalcs=opts.maxCalcs;
maxMatsubara=opts.maxMatsubara;
verbose=opts.verbose;

% Initialising values
nlayers=length(layers);
ninterfaces=nlayers-1+p.cyclic_tunnelling;

% Prefactor             2      e       t    T      a/N^2            /hbar
%                              C       eV  eV      m-2            eV-1 S-1
normalisation_factors = -2*1.60217e-19*p.ts*p.T*p.k_step_size^2/6.582119569e-16;
% Multiplied by E (eV-1) gives C m-2 s-1 = current density

% Building the base matrices with no matsubara frequency dependence
[base_m] = GKTH_hamiltonian(p,layers);
% Matrix for adding frequency dependence
imaginary_identity_m = 1j*eye(4*nlayers,4*nlayers);

% Create matrices for kspace subsampling if you don't already have them
if length(opts.density_grid)~=p.nkpoints
    [density_grid,compute_grid] = GKTH_ksubsample(p,layers);
else
    density_grid=opts.density_grid;
    compute_grid=opts.density_grid;
end
compute_idxs=transpose(find(compute_grid==1));
npointscalc=length(compute_idxs);
%[compute_xs,compute_ys]=ind2sub([p.nkpoints,p.nkpoints],compute_idxs);
random_sampling_max=density_grid(compute_idxs)-1;

% Create the symmetry multiplier
if p.lattice_symmetry == '4mm'
    symmetry_multiplier = 8*ones(p.nkpoints,p.nkpoints);
    for i = 1:p.nkpoints
        symmetry_multiplier(i,i) = 4;
    end
elseif p.lattice_symmetry == 'mm'
    symmetry_multiplier = 4*ones(p.nkpoints,p.nkpoints);
elseif p.lattice_symmetry == 'm'
    symmetry_multiplier = 2*ones(p.nkpoints,p.nkpoints);
else
    symmetry_multiplier = ones(p.nkpoints,p.nkpoints);
end

% Overall multiplier scales the Fs value at each calculated point to
% account for subsampling and symmetry
overall_multiplier = symmetry_multiplier .* density_grid.^2 .* compute_grid;
overall_multiplier_1D = overall_multiplier(compute_idxs)';
%E_matrices_kresolved=zeros(npointscalc,4*nlayers,4*nlayers);
%E_matrices_kresolved_matsum=zeros(npointscalc,4*nlayers,4*nlayers);
%E_matrices=zeros(2,2,ninterfaces);
%E_matrices_new=zeros(2,2,ninterfaces);

px=[0,1;1,0];
py=[0,-1i;1i,0];
pz=[1,0;0,-1];
if include_spin==true
    j_to_include=ones(4*ninterfaces,1);
else
    j_to_include=mod(0:4*ninterfaces-1,4)==0;
end



Expos=4*(1:ninterfaces);
Eypos=Expos-4;
if p.cyclic_tunnelling
    Expos(ninterfaces)=4*(ninterfaces-1);
    Eypos(ninterfaces)=0;
end

% For the first manually defined frequencies
first_freqs=[0,maxMatsubara];
L=length(first_freqs);

%Preallocate array size for speed
ksums=zeros(maxCalcs,ninterfaces*4);
matsubara_freqs=zeros(maxCalcs,1);
weights=zeros(maxCalcs,1);
diffs=zeros(maxCalcs-1,1);
ddiffs=zeros(maxCalcs-1,1);
matsdiffs=zeros(maxCalcs-1,1);
values=zeros(maxCalcs,2);
integrals=zeros(maxCalcs-L,1);

% For verbose only
matsubara_freqs_unsrt=zeros(maxCalcs,1);
E_kresolved=zeros(4,maxCalcs,ninterfaces,npointscalc);
E_kresolved_symm=zeros(4,ninterfaces,p.nkpoints,p.nkpoints);
E_kresolved_final=zeros(4,ninterfaces,2*p.nkpoints,2*p.nkpoints);

% Calculate ksums for first frequencies
matsubara_freqs(1:L)=first_freqs;
for itr=1:L
    ksums(itr,:,:)=calculate_ksum(matsubara_freqs(itr));
    weights(itr)=norm(ksums(itr,:).*j_to_include);
    if verbose==true
        values(itr,:)=[matsubara_freqs(itr) , weights(itr)];
    end
end

% Iterate, dynamically sample new points based on weighting
prev_tol_check=0;
stopped_changing=0;
for itr=L+1:maxCalcs
    % Calculate all the differentials
    diffs(1:itr-2)=(diff(weights(1:itr-1)));
    ddiffs(1:itr-2)=(gradient(diffs(1:itr-2)));
    matsdiffs(1:itr-2)=diff(matsubara_freqs(1:itr-1));
    % Don't try to sample between points already next to each other
    checks=matsdiffs(1:itr-2)>1;
    % Find the new index from the max of the weighting
    [~,mats_idx]=max(abs(checks.*ddiffs(1:itr-2).*matsdiffs(1:itr-2)));
    mats_idx=mats_idx(1); % In case there's two equal values
    new_n=floor((matsubara_freqs(mats_idx)+matsubara_freqs(mats_idx+1))/2);
    % Stick the new n (sorted) into frequency array
    matsubara_freqs(1:itr)=[matsubara_freqs(1:mats_idx);new_n;matsubara_freqs(mats_idx+1:itr-1)];
    % Calculate the sum over k for the new frequency
    new_ksum=calculate_ksum(new_n);
    % Stick this into ksum array
    ksums(1:itr,:)=[ksums(1:mats_idx,:);new_ksum;ksums(mats_idx+1:itr-1,:)];
    weights(1:itr)=[weights(1:mats_idx,:);norm(new_ksum.*j_to_include);weights(mats_idx+1:itr-1)];
    % If verbose, track the integral every step
    if verbose == true
        integrals(itr-L)=trapz(matsubara_freqs(1:itr), weights(1:itr)) + 0.5*weights(1);
        values(itr,:)=[new_n , norm(new_ksum.*j_to_include)];
    end
    % Every 10 iterations check for convergence
    if mod(itr,10)==0
        tol_check = trapz(matsubara_freqs(1:itr), weights(1:itr)) + 0.5*weights(1);
        if abs((tol_check-prev_tol_check)/prev_tol_check) < p.rel_tolerance_Greens
            stopped_changing=stopped_changing+1;
            if stopped_changing==2
                break
            end
        else
            stopped_changing=0;
        end
        prev_tol_check=tol_check;
    end
end
% Finally, integrate over the whole thing to approximate the sum from 0 to
% infinity. You need to add half the first step as trapz cuts this off.
js=zeros(ninterfaces,4);
for interface=1:ninterfaces
    for component=1:4
        idx=sub2ind([4,ninterfaces],component,interface);
        js(interface,component) = trapz(matsubara_freqs, ksums(:,idx) + 0.5*ksums(1,idx));
    end
end

if verbose==true
    [~,id]=sort(matsubara_freqs_unsrt(1:itr));
    E_kresolved=E_kresolved(:,id,:,:);
    for k=1:4
        for j=1:ninterfaces
            for i=1:npointscalc
                E_kresolved_symm(k,j,compute_idxs(i))=trapz(matsubara_freqs(1:itr), E_kresolved(k,1:itr,j,i)) + 0.5*E_kresolved(k,1,j,i);
            end
            E_kresolved_final(k,j,:,:)=GKTH_flipflip(squeeze(E_kresolved_symm(k,j,:,:)),p.use_4mm_symmetry,p.use_kspace_subsampling,density_grid);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%
function [j_txyz] = calculate_ksum(n)
    % Calculate the Matsubara frequency value
    w = (2*n + 1)*pi*p.T;
    ws = imaginary_identity_m*w;
    E_matrices_kresolved=zeros(npointscalc,4*nlayers,4*nlayers);
    % Get the indices for random subsampling
    idx_samps=compute_idxs+round(random_sampling_max.*rand(1,npointscalc))+p.nkpoints*round(random_sampling_max.*rand(1,npointscalc));
    % Go through each k1,k2 point, invert the Hamiltonian to find Fupdown and Fdownup
    for i=1:npointscalc
        idx_samp=idx_samps(i);
        E_matrices_kresolved(i,:,:) = inv(base_m(:,:,idx_samp)+ws) + inv(base_m(:,:,idx_samp)-ws);
    end
    % Finding the new E matrix
    E_matrices=zeros(2,2,ninterfaces);
    j_txyz = zeros(1,ninterfaces*4);
    for i=1:ninterfaces
        for j=1:2
            for k=1:2
                E_matrices(j,k,i) = sum(overall_multiplier_1D.*E_matrices_kresolved(:,Expos(i)+j,Eypos(i)+k)) * normalisation_factors(i);
            end
        end
        j_txyz(4*(i-1)+1) = imag(trace(      E_matrices(:,:,i) ));
        j_txyz(4*(i-1)+2) = imag(trace( px * E_matrices(:,:,i) ));
        j_txyz(4*(i-1)+3) = imag(trace( py * E_matrices(:,:,i) ));
        j_txyz(4*(i-1)+4) = imag(trace( pz * E_matrices(:,:,i) ));
    end
    if verbose==true
        matsubara_freqs_unsrt(itr)=n;
        for n=1:npointscalc
            for i=1:ninterfaces
                for j=1:2
                    for k=1:2
                        E_matrices(j,k,i) = E_matrices_kresolved(n,Expos(i)+j,Eypos(i)+k);
                    end
                end
                E_kresolved(1,itr,i,n) = imag(trace(      E_matrices(:,:,i) ));
                E_kresolved(2,itr,i,n) = imag(trace( px * E_matrices(:,:,i) ));
                E_kresolved(3,itr,i,n) = imag(trace( py * E_matrices(:,:,i) ));
                E_kresolved(4,itr,i,n) = imag(trace( pz * E_matrices(:,:,i) ));
            end
        end
    end
end

end