classdef Global_Params
    %%% Description
    % The Global_Params object contains paramters which apply to a stack of
    % Layers like temperature, applied field and tunnelling paramters, as 
    % well as paramters for controlling the calculation such as function 
    % tolerances and whether to use symmetries or subsampling.
    
    properties
        % Temperature in eV
        T = 1*kB; % this is 1 K
        
        % Global applied field strength
        h = 0;
        % Angle of applied field about y, away from quantization axis z
        theta = 0;
        % Angle of applied field about z, the quantization axis
        theta_ip = 0;
        
        % Tunnelling parameters. Should be length nlayers-1.
        ts=zeros(100,1)+0.0005;
        % Whether the last layer tunnels to the first layer
        cyclic_tunnelling=false;
        
        % lattice parameter
        a = 3.905e-10;
        
        % direction of broken inversion symmetry. x=1,y=2,z=3
        interface_normal=3;
        
        % How big to make the area in k-space. 1 goes to 1st Brillouin
        % zone, 2 goes to the 3rd Brillouin zone, 3 goes to the 6th
        % Brillouin zone.
        BZ_multiplier=1;
        
        
        % Number of kpoints: must be multiple of 8 for meshing to work with mm symmetry
        % Minimum 80 for testing, 160 or greater for real use
        nkpoints = 160;
        
        % Symmetries and tricks for code accelerating.
        %use_m_symmetry=false;
        m_symmetry_line=0; % angle of the line of symmetry
        %use_mm_symmetry=true;
        %use_4mm_symmetry=true;
        
        % Would be better, not currently used
        lattice_symmetry (1,1) string {mustBeMember(lattice_symmetry,{'4mm','mm','m','none'})} ='4mm';
        
        % Not useful if using radials
        use_kspace_subsampling=true;
        subsampling_point_fraction=0.25;
        
        % I don't think you can disable this now
        use_matsubara_subsampling=true;
        
        % For the radial calculations
        nradials (1,1) {mustBeInteger} = 50;
        ntest (1,1) {mustBeInteger} = 2000;
        nfinal (1,1) {mustBeInteger} = 200;
        
        
        % Tolerances for Greens function iterations etc.
        
        % Tolerances for the Greens function sum over matsubara frequencies
        % These get multiplied by T in (eV) to account for slower
        % convergence at low temperature.
        abs_tolerance_Greens=1e-99; %This needs to be low so the hc finder works. Might as well not exist.
        rel_tolerance_Greens=1e-5;
        
        abs_tolerance_hc=1e-8;
        
        % Absolute tolerance used in fzero self consistency
        abs_tolerance_self_consistency_1S=1e-6;
        
        % Relative tolerance used in iterative self consistency
        abs_tolerance_self_consistency_iterate=1e-7;
        rel_tolerance_self_consistency_iterate=1e-6;
        
    end
    
    properties (Dependent)
        % k1 and k2 are nxn arrays, defined by obj.nkpoints, and contains 
        % k1 or k2 values at each point. The diagonals are shifted slightly 
        % to avoid divide-by-zeros when k1=k2.
        k1
        k2
        k_step_size
    end
    
    methods
        
        function obj = set.nkpoints(obj,value)
            value=floor(value);
            obj.nkpoints = value + 15 - mod(value-1,16);
        end
        
        function k1 = get.k1(obj)
            scale=obj.BZ_multiplier*pi/obj.a;
            if obj.lattice_symmetry == "mm" || obj.lattice_symmetry == "4mm"
                ks=scale*linspace((1/(2*obj.nkpoints)),(1-1/(2*obj.nkpoints)),obj.nkpoints);
                [k1,~]=meshgrid(ks,ks);
%                 noisify=(rand(obj.nkpoints)-0.5)*scale/(obj.nkpoints);
%                 noisify=noisify.*~eye(obj.nkpoints);
%                 k1=k1+noisify;
                for i=1:obj.nkpoints
                    k1(i,i)= k1(i,i)+(2*mod(i,2)-1)*scale/(obj.nkpoints*10); %shift by a 10th of a step size with alternating sign to avoid divide by zeros
                end
            else
                ks=scale*linspace((-1+1/(obj.nkpoints)),(1-1/(obj.nkpoints)),obj.nkpoints);
                [k1,~]=meshgrid(ks,ks);
                for i=1:obj.nkpoints
                    k1(i,i)= k1(i,i)+(2*mod(i,2)-1)*scale/(obj.nkpoints*10); %shift by a 10th of a step size with alternating sign
                    k1(i,1+obj.nkpoints-i) = k1(i,1+obj.nkpoints-i)+(2*mod(i,2)-1)*scale/(obj.nkpoints*10); %shift by a 5th of a step size with alternating sign
                end
            end
        end
        
        function k2 = get.k2(obj)
            scale=obj.BZ_multiplier*pi/obj.a;
            if obj.lattice_symmetry == "mm" || obj.lattice_symmetry == "4mm"
                ks=scale*linspace((1/(2*obj.nkpoints)),(1-1/(2*obj.nkpoints)),obj.nkpoints);
                [~,k2]=meshgrid(ks,ks);
%                 noisify=(rand(obj.nkpoints)-0.5)*scale/(obj.nkpoints);
%                 noisify=noisify.*~eye(obj.nkpoints);
%                 k2=k2+noisify;
                for i=1:obj.nkpoints
                    k2(i,i)= k2(i,i)-(2*mod(i,2)-1)*scale/(obj.nkpoints*10); %shift by a 10th of a step size with alternating sign
                end
            else
                ks=scale*linspace((-1+1/(obj.nkpoints)),(1-1/(obj.nkpoints)),obj.nkpoints);
                [~,k2]=meshgrid(ks,ks);
                for i=1:obj.nkpoints
                    k2(i,i)= k2(i,i)-(2*mod(i,2)-1)*scale/(obj.nkpoints*10); %shift by a 10th of a step size with alternating sign
                    k2(i,1+obj.nkpoints-i) = k2(i,1+obj.nkpoints-i)+(2*mod(i,2)-1)*scale/(obj.nkpoints*10); %shift by a 10th of a step size with alternating sign
                end
            end
        end
        
        function k_step_size = get.k_step_size(obj)
            if obj.lattice_symmetry == "mm" || obj.lattice_symmetry == "4mm"
                k_step_size = obj.BZ_multiplier*pi/(obj.a*obj.nkpoints);
            else
                k_step_size = 2*obj.BZ_multiplier*pi/(obj.a*obj.nkpoints);
            end
        end
    end
end