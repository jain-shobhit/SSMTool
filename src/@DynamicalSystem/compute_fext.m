function fext = compute_fext(obj,t)
% COMPUTE_FEXT We compute the external force at a particular time t 
% in a second-order mechanical system. 

assert(obj.order == 2, ' fext can only be computed for second-order systems')

if isempty(obj.fext)
    fext = sparse(obj.n,1);
else
    if obj.Options.HarmonicForce
        fext = obj.fext.epsilon * real(obj.fext.coeffs * exp(1i * obj.fext.kappas * obj.Omega * t));
        if obj.Options.BaseExcitation
            fext = fext*(obj.Omega)^2; 
        end
    else
        fext = obj.fext(t);
    end    
end
