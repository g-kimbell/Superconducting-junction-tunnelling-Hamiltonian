function [eigenvalues] = GKTH_find_spectrum(p,layers)

%%% Description
% Finds the eigenvalues of a layer stack, which correspond to the band
% energies
%
%%% Inputs
%   p      :    Global Parameter object defining the stack and 
%               calculation paramters.
%   layers :    An array of Layer objects defining the stack
%
%%% Outputs
% eigenvalues:  n x n x 4*layers matrix of k-resolved eigenvalues, with no
%               particular order along the 3rd axis

% Set all Deltas to zero and get Hamiltonian
nlayers=length(layers);
temp=num2cell(zeros(nlayers,1));
[layers.Delta_0]=temp{:};
[Hs] = GKTH_hamiltonian(p,layers);

% This just leaves them sorted in the way matlab guesses
eigenvalues=zeros(p.nkpoints,p.nkpoints,4*nlayers);
for i=1:p.nkpoints
    for j=1:p.nkpoints
        eigenvalues(i,j,:)=eig(Hs(:,:,i,j));
    end
end