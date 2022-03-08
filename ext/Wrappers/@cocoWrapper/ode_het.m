function y = ode_het(obj, t, z, p, data)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - omega

if isempty(data)
    obj.system.Omega        = p(1:end-1,:);
    obj.system.fext.epsilon = p(end,:);
    y = obj.system.odefun(t, z);
else
    n = obj.system.n;
    nt = size(z,2);
    x  = z(1:n,:); 
    xd = z(n+1:2*n,:);
    om = p(1:end-1,:);
    ep = p(end,:);

    % linear part
    y1 = obj.system.BinvA*z;

    % nonlinear part
    fnl = data.fnl; % fnl.coeffs and fnl.ind
    numNonlinearTerms = size(fnl.coeffs,2);
    y2 = 0;
    for i=1:numNonlinearTerms
        coeff = repmat(fnl.coeffs(:,i),[1, nt]);
        ind = fnl.ind(i,:);
        % find nonzero exponents
        expind = find(ind);
        s = 1;
        for j=1:numel(expind)
            s = s.*z(expind(j),:).^ind(expind(j));
        end
        s = repmat(s, [n, 1]);
        y2 = y2+coeff.*s;
    end
    
    % external forcing
    assert(~isempty(obj.system.fext), 'no external forcing');
    fext_coeffs = repmat(obj.system.fext.coeffs(:,1), [1, nt]);
    if data.isbaseForce
        fext_harm = repmat(ep.*om.^2.*cos(obj.system.fext.kappas(1)*om.*t), [n, 1]);
    else
        fext_harm = repmat(ep.*cos(obj.system.fext.kappas(1)*om.*t), [n, 1]);
    end
    y3 = 2*fext_coeffs.*fext_harm;
    
    y = y1 + [zeros(n,nt); obj.system.M\(-y2+y3)];
end
end