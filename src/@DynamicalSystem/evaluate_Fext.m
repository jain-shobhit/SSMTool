function Fext = evaluate_Fext(obj,t)
% EVALUATE_FEXT This function evaluates the external forcing at a given
% time t for a dynamical system object in the first-order form 

switch obj.order
    case 1
        if isempty(obj.Fext)
            Fext = sparse(obj.N,1);
        else
            Fext = obj.Fext.epsilon * real(obj.Fext.coeffs * exp(1i * obj.Fext.kappas * obj.Omega * t));
        end
        
    case 2
        Fext = [obj.compute_fext(t);
                sparse(obj.n,1)];
end