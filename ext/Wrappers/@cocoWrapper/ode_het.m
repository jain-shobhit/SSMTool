function y = ode_het(obj, t, z, p)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - omega

om = p(1:end-1,:);
ep = p(end,:);

Omega = obj.system.Omega;
obj.system.Omega = om;
switch obj.system.order
    case 1
        % record original values
        epsilon = obj.system.Fext.epsilon;
        % assign vector values
        obj.system.Fext.epsilon = ep;
        % evaluate RHS
        y = obj.system.odefun(t,z);
        % assign back orignal values
        obj.system.Fext.epsilon = epsilon;
    case 2
    n = obj.system.n;
    x  = z(1:n,:); 
    xd = z(n+1:2*n,:);
        % record original values
        epsilon = obj.system.fext.epsilon;
        % assign vector values
        obj.system.fext.epsilon = ep;
        % evaluate RHS
        y = obj.system.odefun(t,[x;xd]);
        % assign back orignal values
        obj.system.fext.epsilon = epsilon;
    otherwise
        error('Dynamical system order must be 1 or 2')
end
obj.system.Omega = Omega;
