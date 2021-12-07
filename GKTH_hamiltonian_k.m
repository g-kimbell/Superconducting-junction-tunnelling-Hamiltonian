function[m] = GKTH_hamiltonian_k(p,k1s,k2s,layers)

%%% Description
% Creates a (4*nlayers)*(4*nlayers) Hamiltonian matrix for every pair of
% k-points defined by the matrices k1s and k2s
%
%%% Inputs
%   p           :   Global Parameter object defining the stack and 
%                   calculation paramters
%   k1s         :   A 2D array of kx points
%   k2s         :   A 2D array of ky points
%   layers      :   An array of Layer objects, defining the junction 
%                   structure.
%
%%% Outputs
%   m           :   The Hamiltonian matrix defined at every k-point

nks=length(k1s);
nlayers=length(layers);
%%% Base matrix %%%
m = zeros([(4*nlayers) (4*nlayers) size(k1s)]);

function m = rotate_y(m,theta)
    Ry = [cos(theta/2) sin(theta/2) ; -sin(theta/2) cos(theta/2)];
    Ryt = ctranspose(Ry);
    m = Ry*m*Ryt;
end
function m = rotate_z(m,theta_ip)
    Rz = [exp(1i*theta_ip/2) 0 ; 0 exp(-1i*theta_ip/2)];
    Rzt = ctranspose(Rz);
    m = Rz*m*Rzt;
end

%%% Fill in each layer along the diagonal %%%
for i=1:nlayers
    L=layers(i);
    idx=(i-1)*4;
    %%% Electronic spectrum
    xis = zeros([2,2,size(k1s)]);
    xis(1,1,:,:) = -GKTH_spectrum_k(p,L,k1s,k2s);
    xis(2,2,:,:) = -GKTH_spectrum_k(p,L,k1s,k2s);
    %%% Magnetism
    h_layer  = rotate_z(rotate_y([L.h 0 ; 0 -L.h-L.dE ],-L.theta),-L.theta_ip);
    h_global = rotate_z(rotate_y([p.h 0 ; 0   -p.h    ],-p.theta),-p.theta_ip);
    %%% SOC
    if L.alpha ~=0
        SOC=zeros([2,2,size(k1s)]);
        switch p.interface_normal
            case 1
                SOC(1,1,:,:) = -k1s;
                SOC(1,2,:,:) = -1j*k2s;
                SOC(2,1,:,:) = 1j*k2s;
                SOC(2,2,:,:) = k1s;
            case 2
                SOC(1,1,:,:) = k2s;
                SOC(1,2,:,:) = -k1s;
                SOC(2,1,:,:) = -k1s;
                SOC(2,2,:,:) = -k2s;
            case 3
                SOC(1,2,:,:) = 1j*k1s + k2s;
                SOC(2,1,:,:) = -1j*k1s + k2s;
        end
        SOC = - L.alpha * p.a/pi * SOC;
    else
        SOC=0;
    end
    m(1+idx:2+idx,1+idx:2+idx,:,:) = xis + h_layer + h_global + SOC;
    m(3+idx:4+idx,3+idx:4+idx,:,:) =  - conj(xis + h_layer + h_global + SOC);
    
    %%% Superconductivity
    m(1+idx,4+idx,:,:) =  exp(1j*L.phi)*GKTH_Delta_k(p,L,k1s,k2s);
    m(2+idx,3+idx,:,:) =  -exp(1j*L.phi)*GKTH_Delta_k(p,L,k1s,k2s);
    m(3+idx,2+idx,:,:) =  -exp(-1j*L.phi)*GKTH_Delta_k(p,L,k1s,k2s);
    m(4+idx,1+idx,:,:) =  exp(-1j*L.phi)*GKTH_Delta_k(p,L,k1s,k2s);
end

%%% Tunnelling between layers %%%
signs=[-1,-1,1,1];
for i=1:(nlayers-1)
    idx=(i-1)*4;
    for j=1:4
    	m(j+idx,j+4+idx,:,:)=p.ts(i)*signs(j);
        m(j+4+idx,j+idx,:,:)=p.ts(i)*signs(j);
    end
end
if p.cyclic_tunnelling==true
    shift=4*(nlayers-1);
    for j=1:4
        m(j,j+shift,:,:)=p.ts(i)*signs(j);
        m(j+shift,j,:,:)=p.ts(i)*signs(j);
    end
end

end