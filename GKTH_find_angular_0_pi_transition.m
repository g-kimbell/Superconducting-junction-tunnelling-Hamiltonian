function[theta_zero] = GKTH_find_angular_0_pi_transition(p,layers,interval,opts)

%%% Description
% Function to find 0-pi transitions as the magnetisation vector of one or
% more layers of a junction is rotated.
%
%%% Inputs
%   p      :    Global Parameter object defining the stack and 
%               calculation paramters.
%   layers :    An array of Layer objects defining the stack
%  interval:    The interval in theta over which to check for 0-pi 
%               transitions
%   opts   :    Optional keyword arguements
%      tolerance: the absolute tolerance for finding the angle of the 0-pi
%                 transition
%      layer_to_rotate: array of integers defining which layers are rotated
%      layer_to_phase: which layer to change the phase of
%      in_plane: bool, whether to rotate the magnetisation in-plane (true)
%                or out-of-plane (false)
%      type: string defining the type of rotation, either 'sym' for every
%            layer rotating the same way, 'asym' for layers rotating in
%            alternating directions, or 'helix' for gradual rotation of
%            layers.
%          
%%% Outputs
% theta_zero: the theta value at which there is a 0-pi transition
%

arguments
    p
    layers
    interval
    opts.tolerance = 1e-4;
    opts.layer_to_rotate=[1,3];
    opts.layer_to_phase=3;
    opts.in_plane=true;
    opts.type (1,1) string {mustBeMember(opts.type,{'sym','asym','helix'})} = "asym";
end

function rotate(theta)
    for i=1:length(opts.layer_to_rotate)
        if opts.type=="sym"
            multiply=1;
        elseif opts.type=="asym"
            multiply=mod(i,2)*2-1;
        elseif opts.type=="helix"
            multiply=i;
        end
        if opts.in_plane==true
            layers(opts.layer_to_rotate(i)).theta_ip=multiply*theta;
        else
            layers(opts.layer_to_rotate(i)).theta=multiply*theta;
        end
    end
end

min_theta=interval(1);
max_theta=interval(2);
layers(opts.layer_to_phase).phi=pi/2;

%First check theres a sign change
rotate(min_theta)
[js] = GKTH_Greens_current_radial(p,layers);
min_j=js(1);
rotate(max_theta);
[js] = GKTH_Greens_current_radial(p,layers);
max_j=js(1);
if min_j*max_j < 0 %there's a sign change
    dtheta=abs(max_theta-min_theta);
    while dtheta>opts.tolerance
        % This approach is sometimes fast but sometimes not robust
        %theta_guess=interp1([min_j,max_j],[min_theta,max_theta],0);
        theta_guess=(max_theta+min_theta)/2;
        rotate(theta_guess);
        [js] = GKTH_Greens_current_radial(p,layers);
        if js(1)*min_j<0 % opposite sign
            max_j=js(1);
            max_theta=theta_guess;
        else
            min_j=js(1);
            min_theta=theta_guess;
        end
        dtheta=abs(max_theta-min_theta);
    end
    theta_zero=theta_guess;
else
    theta_zero=NaN;
end

end