function f = odefun(obj,t,z)
% ODEFUN this function computes the right-hand-side of the ODE 
% described by the DynamicalSystem object (obj)

switch obj.order
    case 1
        f = obj.B\(obj.A * z + obj.evaluate_Fnl(z) + obj.Fext(t));
    case 2
        x = z(1:obj.n);
        xd = z(obj.n+1:obj.N);
        f = obj.BinvA * z + ...
            [sparse(obj.n,1); 
            obj.M\(-obj.compute_fnl(x,xd) + obj.compute_fext(t))];
end
end