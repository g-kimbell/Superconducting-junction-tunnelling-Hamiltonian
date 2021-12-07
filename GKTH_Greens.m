function [Fs_sums,matsubara_freqs,ksums,integrals,values,F_kresolved_final] = GKTH_Greens(p,layers,opts)

%%% Description
% GKTH_Greens inverts the Hamiltonian and takes a sum over k and 
% Matsubara frequencies, returning the Fs sum for each layer.
%
%%% Inputs:
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
%   opts        :   Optional keyword arguements
%        density_grid : density_grid from GKTH_ksubsample if you know the 
%                       Fermi surface won't change during repeated Greens 
%                       calls. If unsure don't pass it and it will
%                       calculate a new one each time.
%        compute_grid : compute_grid from GKTH_ksubsample if you know the 
%                       Fermi surface won't change during repeated Greens
%                       calls, if unsure leave blank.
%        maxCalcs     : Maximum number of Matsubara calculations
%        maxMatsubara : The maximum Matsubara frequency used. May need
%                       increasing if you are looking at very low 
%                       temperatures ( < 0.1 K )
%        verbose      : Whether or not to store and return the k-resolved 
%                       anomalous Green function and the integrated values 
%                       vs Matsubara frequency
%
%%% Outputs:
%
%   Fs_sums        :  Array containing the Fs sum over k1, k2 and 
%                     Matsubara frequencies for each layer in the stack.
% matsubara_freqs  :  The Matsubara frequencies calculated
%   ksums          :  The total sum over k for each frequency calculated
% integrals        :  (verbose=true only) The k-sum integrated over
%                     Matsubara frequency vs matsubara frequency
%   values         :  (verbose=true only) The values used in the weighting 
%                     function
% F_kresolved_final:  (verbose=true only) The final anomalous Green 
%                     function integrated over Matsubara frequency resolved
%                     in k-space
%

arguments
    p
    layers
    opts.density_grid = NaN;
    opts.compute_grid = NaN;
    opts.maxCalcs=500;
    opts.maxMatsubara=1e7; % before 1e6
    opts.verbose=false;
end
maxCalcs=opts.maxCalcs;
maxMatsubara=opts.maxMatsubara;
verbose=opts.verbose;

% Initialising values
nlayers=length(layers);
D_factors=zeros(p.nkpoints,p.nkpoints,nlayers);
for i=1:nlayers
    D_factors(:,:,i)=1./(GKTH_Delta(p,layers(i),1));
end
D_factors(isnan(D_factors))=0;
normalisation_factor=p.k_step_size^2/(2*pi/p.a)^2;
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
overall_multiplier = overall_multiplier(compute_idxs);

% For the first manually defined frequencies
first_freqs=[0,maxMatsubara];
L=length(first_freqs);

%Preallocate array size for speed
ksums=zeros(maxCalcs,nlayers);
matsubara_freqs=zeros(maxCalcs,1);
weights=zeros(maxCalcs,1);
diffs=zeros(maxCalcs-1,1);
ddiffs=zeros(maxCalcs-1,1);
matsdiffs=zeros(maxCalcs-1,1);
values=zeros(maxCalcs,1+nlayers);
integrals=zeros(maxCalcs-L,1);

% For verbose only
matsubara_freqs_unsrt=zeros(maxCalcs,1);
F_kresolved=zeros(maxCalcs,nlayers,npointscalc);
F_kresolved_symm=zeros(nlayers,p.nkpoints,p.nkpoints);
F_kresolved_final=zeros(nlayers,2*p.nkpoints,2*p.nkpoints);

% Calculate ksums for first frequencies
matsubara_freqs(1:L)=first_freqs;
for itr=1:L
    ksums(itr,:)=calculate_ksum(matsubara_freqs(itr));
    if verbose==true
        values(itr,:)=[matsubara_freqs(itr) , abs(ksums(itr,:))];
    end
end
weights(1:L)=sum(abs(ksums(1:L,:)),2);
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
    weights(1:itr)=[weights(1:mats_idx,:);sum(abs(new_ksum));weights(mats_idx+1:itr-1)];
    
    % If verbose, track the integral every step
    if verbose == true
        integrals(itr-L)=trapz(matsubara_freqs(1:itr), weights(1:itr)) + 0.5*weights(1);
        values(itr,:)=[new_n , abs(new_ksum)];
    end
    % Every 10 iterations check for convergence
    if mod(itr,10)==0
        tol_check = trapz(matsubara_freqs(1:itr), weights(1:itr)) + 0.5*weights(1);
        if abs((tol_check-prev_tol_check)/prev_tol_check) < p.rel_tolerance_Greens
            stopped_changing=stopped_changing+1;
            if stopped_changing==3
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
Fs_sums=zeros(nlayers,1);
for layer=1:nlayers
    Fs_sums(layer) = trapz(matsubara_freqs(1:itr), ksums(1:itr,layer)) + 0.5*ksums(1,layer);
end

if verbose==true
    [~,id]=sort(matsubara_freqs_unsrt(1:itr));
    F_kresolved=F_kresolved(id,:,:);
    for j=1:nlayers
        for i=1:npointscalc
            F_kresolved_symm(j,compute_idxs(i))=trapz(matsubara_freqs(1:itr), F_kresolved(1:itr,j,i)) + 0.5*F_kresolved(1,j,i);
        end
        F_kresolved_final(j,:,:)=GKTH_flipflip(squeeze(F_kresolved_symm(j,:,:)),p.use_4mm_symmetry,p.use_kspace_subsampling,density_grid);
        
    end
    subplot(1,2,1)
    pcolor(real(squeeze(F_kresolved_final(1,:,:))))
    subplot(1,2,2)
    plot(matsubara_freqs,abs(ksums),'o-')
end

%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%
function [Fs_ksum] = calculate_ksum(n)
    Fs_ksum=zeros(1,nlayers);
    w = (2*n + 1)*pi*p.T;
    ws = imaginary_identity_m*w;
    Fupdowns=zeros(nlayers,npointscalc);
    Fdownups=zeros(nlayers,npointscalc);
    % Get the indices for random subsampling
    idx_samps=compute_idxs+round(random_sampling_max.*rand(1,npointscalc))+p.nkpoints*round(random_sampling_max.*rand(1,npointscalc));
    % Go through each k1, k2 point, invert the Hamiltonian to find Fupdown and Fdownup
    for i=1:npointscalc
        idx_samp=idx_samps(i);
        m_inv = inv(base_m(:,:,idx_samp)+ws);
        for j=0:nlayers-1
            Fupdowns(j+1,i)=m_inv(1+j*4,4+j*4)*D_factors(idx_samp+j*p.nkpoints^2);
            Fdownups(j+1,i)=m_inv(2+j*4,3+j*4)*D_factors(idx_samp+j*p.nkpoints^2);
        end 
    end
    % For sums multiply by 2 for frequency symmetry, divide by 2 for singlet Fs=1/2 (Fupdown-Fdownup)
    for i=1:nlayers
        Fs_ksum(i)= sum(overall_multiplier.*(Fupdowns(i,:) - Fdownups(i,:))) * normalisation_factor;
        if verbose==true
            F_kresolved(itr,i,:)=Fupdowns(i,:) - Fdownups(i,:);
            matsubara_freqs_unsrt(itr)=n;
        end
    end
    
end

end