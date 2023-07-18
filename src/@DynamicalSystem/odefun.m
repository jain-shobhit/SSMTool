function f = odefun(obj,t,z)
% ODEFUN this function computes the right-hand-side of the ODE 
% described by the DynamicalSystem object (obj)
        f = obj.B\(obj.A * z + obj.evaluate_Fnl(z) + obj.evaluate_Fext(t,z));
end