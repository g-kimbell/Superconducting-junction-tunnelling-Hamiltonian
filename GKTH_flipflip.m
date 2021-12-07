function [new_m] = GKTH_flipflip(m,use_4mm_symmetry,use_kspace_subsampling,density_grid)

%%% Description:
% GKTH_flipflip constructs a full matrix out of a sub-calculated matrix
% from the symmetry used. E.g. the new Fs sum is often calculated only on
% 1/8 of the Brillouin zone if 4mm symmetry is used. Flip flip
% reconstrcucts the rest of the Brillouin zone from this matrix.
%
%%% Inputs:
%   m                   :   The input matrix to reconstruct
%   use_4mm_symmetry    :   Whether 4mm symmetry was used to construct the
%                           matrix m.
% use_kspace_subsampling:   Whether grid-base k-space subsampling was used
%                           in the calculation
% density_grid          :   The density grid for k-space subsampling if
%                           used
%
%%% Outputs:
%   new_m               :   The input matrix after reconstruction
%

arguments
    m
    use_4mm_symmetry=false
    use_kspace_subsampling=false
    density_grid=ones(length(m))
end

len=length(m);

if use_kspace_subsampling == true
    for i = 1:len
        for j = i:len
            d=density_grid(i,j);
            m(i,j)=m(floor((i-1)/d)*d+1,floor((j-1)/d)*d+1);
        end
    end
end

if use_4mm_symmetry == true
    for i = 1:len
        for j = i:len
            m(j,i)=m(i,j);
        end
    end
end

new_m=zeros([2*len,2*len]);
new_m(len+1:2*len,len+1:2*len)=m;
new_m(1:len,len+1:2*len)=flip(m,1);
new_m(len+1:2*len,1:len)=flip(m,2);
new_m(1:len,1:len)=flip(flip(m,1),2);