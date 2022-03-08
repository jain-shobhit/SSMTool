function Fext = evaluate_Fext(obj,t)
% EVALUATE_FEXT This function evaluates the external forcing at a given
% time t for a dynamical system object in the first-order form 

switch obj.order
    case 1
        if isempty(obj.Fext)
            Fext = sparse(obj.N,1);
        else
            if obj.Options.HarmonicForce
                Fext = obj.Fext.epsilon * real(obj.Fext.coeffs * exp(1i * obj.Fext.kappas * obj.Omega * t));
                if obj.Options.BaseExcitation
                    Fext = Fext*(obj.Omega)^2;
                end
            else
                Fext = obj.Fext(t);
            end
        end
        
    case 2
        Fext = [obj.compute_fext(t);
                sparse(obj.n,1)];
end
