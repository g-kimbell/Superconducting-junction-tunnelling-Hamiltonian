classdef Layer
    %%% Description
    % A Layer object contains all of its properties which are used to
    % define its electronic spectrum, gap, and its contribution to the
    % overall Hamiltonian.
    
    properties  
        %%%%% Superconductivity
        % Delta_0
        Delta_0 (1,1) {mustBeNumeric,mustBeReal} = 0;
        % Symmetry of the superconductor: s=swave d=dwave or n=normal
        symmetry = "s";
        % Superconducting coupling strength
        lambda (1,1) {mustBeNumeric,mustBeReal} = 0;
        % Superconducting phase
        phi (1,1) {mustBeNumeric,mustBeReal} = 0;

        %%%%%% Ferromagnetism
        % Symmetric exchange field in eV acting on both up and down
        h (1,1) {mustBeNumeric,mustBeReal} = 0;
        % Asymmetric exchange shift of only one band in eV
        dE (1,1) {mustBeNumeric,mustBeReal} = 0;
        % Angle of exchange field about y, away from quantization axis z
        theta (1,1) {mustBeNumeric,mustBeReal} = 0;
        % Angle of exchange field about z, the quantization axis
        theta_ip (1,1) {mustBeNumeric,mustBeReal} = 0;
        
        %%%%%% SOC
        alpha = 0;

        %%%%%% Electronic spectrum
        % Density of states at Fermi surface, just normalised to 1 for now
        N0 (1,1) {mustBeNumeric,mustBeReal} = 1;
        % Nearest neighbour hopping parameter in eV
        tNN (1,1) {mustBeNumeric,mustBeReal} = -0.7823;
        % Next-nearest neighbour hopping parameter in eV
        tNNN (1,1) {mustBeNumeric,mustBeReal} = -0.0740;
        % Chemical potential in eV
        mu (1,1) {mustBeNumeric,mustBeReal} = 0.06525;
        % Type of dispersion
        dispersion_type = "tb";
    end
    
    methods
        
        %%% Superconducting gap
        function obj = Layer(n)
            if nargin ~= 0
                obj(n)=obj;
            end
        end
        
        %%% Superconducting gap
        function Ds = Ds(obj,p)
            Ds=GKTH_Delta(p,obj,obj.Delta_0);
        end
        
        %%% Electronic spectrum
        function xis = xis(obj,p)
            if strcmp(obj.dispersion_type,'tb')
                xis = obj.tNN*(cos(p.k1*p.a)+cos(p.k2*p.a)) + obj.tNNN*cos(p.k1*p.a).*cos(p.k2*p.a) + obj.mu;
            elseif strcmp(obj.dispersion_type,'para')
                hbar=6.582119569e-16;
                me=9.10938356e-31;
                xis = hbar^2*((p.k1*p.a).^2+(p.k2*p.a).^2)/(2*me)+obj.mu;
            else
                error("Spectrum type must be 'tb' for tight binding or 'para' for parabolic.")
            end
        end
        
    end
end