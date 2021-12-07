function [Fs_sums,F_kresolved_matsum,k1s,k2s] = GKTH_Greens_radial(p,layers,opts)

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
%        maxCalcs     : Maximum number of Matsubara calculations
%        maxMatsubara : The maximum Matsubara frequency used. May need
%                       increasing if you are looking at very low 
%                       temperatures ( < 0.1 K )
%        verbose      : Whether or not to store and return the k-resolved 
%                       anomalous Green function and the integrated values 
%                       vs Matsubara frequency
%       just_use_layer: Define which layer(s) to use for the radial k-space
%                       sub-sampling
%
%%% Outputs:
%
%   Fs_sums        :  Array containing the Fs sum over k1, k2 and 
%                     Matsubara frequencies for each layer in the stack.
% F_kresolved_matsum: (verbose=true only) The final anomalous Green 
%                     function integrated over Matsubara frequency resolved
%                     in k-space
%  k1s             :  The kx points for the k-resolved anomalous Green
%                     function
%  k1s             :  The ky points for the k-resolved anomalous Green
%                     function

arguments
    p
    layers
    opts.maxCalcs=500;
    opts.maxMatsubara=1e7; % before 1e6
    opts.verbose=false;
    opts.just_use_layer=1;
end
maxCalcs=opts.maxCalcs;
maxMatsubara=opts.maxMatsubara;
verbose=opts.verbose;

% Initialising values
nlayers=length(layers);

[k1s,k2s,~,~,area_factor] = GKTH_find_radial_ks(p,layers,just_use_layer=opts.just_use_layer);

[nrs,nangles]=size(k1s);

D_factors=zeros(nlayers,nrs,nangles);
for i=1:nlayers
    L=layers(i);
    L.Delta_0=1;
    D_factors(i,:,:)=1./GKTH_Delta_k(p,L,k1s,k2s);
end
D_factors(isinf(D_factors)|isnan(D_factors))=0;

% Building the base matrices with no matsubara frequency dependence
%[base_m] = GKTH_hamiltonian_k(p,k1s,k2s,layers);
[base_m] = GKTH_hamiltonian_k(p,k1s,k2s,layers);

% Matrix for adding frequency dependence
imaginary_identity_m = 1j*eye(4*nlayers,4*nlayers);

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
values=zeros(maxCalcs,nlayers+1);
integrals=zeros(maxCalcs-L,1);

% For verbose only
matsubara_freqs_unsrt=zeros(maxCalcs,1);
F_kresolved=zeros(maxCalcs,nlayers,nrs,nangles);
F_kresolved_matsum=zeros(nlayers,nrs,nangles);

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
Fs_sums=zeros(nlayers,1);
for layer=1:nlayers
    Fs_sums(layer) = trapz(matsubara_freqs(1:itr), ksums(1:itr,layer)) + 0.5*ksums(1,layer);
end

if verbose==true
   [~,id]=sort(matsubara_freqs_unsrt(1:itr));
   F_kresolved=F_kresolved(id,:,:,:);

   for j=1:nlayers
       for i1=1:nrs
           for i2=1:nangles
                F_kresolved_matsum(j,i1,i2)=trapz(matsubara_freqs(1:itr), F_kresolved(1:itr,j,i1,i2)) + 0.5*F_kresolved(1,j,i1,i2);
           end
       end
   end
end


%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%
function [Fs_ksum] = calculate_ksum(n)
    % Calculate the Matsubara frequency value
    Fs_ksum=zeros(1,nlayers);
    w = (2*n + 1)*pi*p.T;
    ws = imaginary_identity_m*w;
    Fupdowns=zeros(nlayers,nrs,nangles);
    Fdownups=zeros(nlayers,nrs,nangles);
    % Go through each k1,k2 point, invert the Hamiltonian to find Fupdown and Fdownup
    for i1=1:nrs
        for i2=1:nangles
            m_inv = inv(base_m(:,:,i1,i2)+ws);
            for j=1:nlayers
                Fupdowns(j,i1,i2)=m_inv(1+(j-1)*4,4+(j-1)*4);
                Fdownups(j,i1,i2)=m_inv(2+(j-1)*4,3+(j-1)*4);
            end 
        end
    end
    % For sums multiply by 2 for frequency symmetry, divide by 2 for singlet Fs=1/2 (Fupdown-Fdownup)
    for i=1:nlayers
        Fs_ksum(i)= sum(area_factor .* squeeze(D_factors(i,:,:) .* (Fupdowns(i,:,:) - Fdownups(i,:,:))),'all');
        if verbose==true
            F_kresolved(itr,i,:,:) = squeeze(D_factors(i,:,:).*(Fupdowns(i,:,:) - Fdownups(i,:,:)));
        end
    end
    if verbose==true
       matsubara_freqs_unsrt(itr)=n;
    end
end

end