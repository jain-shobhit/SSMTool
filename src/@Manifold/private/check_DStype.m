function [DStype] = check_DStype(obj)
% Check if A,B or the nonlinearity feature complex coefficients - if yes,
% then there are no symmetries that can be exploited for computation.


if ~isreal(obj.System.A) 
    obj.System.Options.DStype  = 'complex';
    
elseif ~isreal(obj.System.B)
    obj.System.Options.DStype  = 'complex';
    
elseif ~checkF(obj)  
    obj.System.Options.DStype  = 'complex';
    
elseif obj.dimManifold == 1
    % System is not complex but computational routine works exactly the
    % same for an overdamped master mode as for complex DS
    obj.System.Options.DStype  = 'complex';
    
elseif obj.System.Options.ChooseComplexComp
    obj.System.Options.DStype  = 'complex';
elseif ~obj.System.Options.ChooseComplexComp
    obj.System.Options.DStype  = 'real';
end

DStype = obj.System.Options.DStype;
end


function [indicator] = checkF(obj)


d = length(obj.System.F) ;
for j = 1:d
    if ~isreal(obj.System.F(j).coeffs)
        indicator = false;
        return
    end
end

indicator = true;
end
