function y = ode_aut(obj, z, p, data)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - (epsilon, alpha) where epsilon is a parameter will
% be applied to damping matrix, and alpha is the parameter will be applied
% to nonlinearity. We assume there is only one conserved quantity in the
% system without damping. If the system has other symmetries, this
% implementation will fail.

if isempty(data)
    error('a data structure with fnl field should be given');
else
    n = obj.system.n;
    nt = size(z,2);
    x  = z(1:n,:); 
    xd = z(n+1:2*n,:);

    % linear part
    y1x = xd;
    y1d = obj.system.K*x+repmat(p(1,:),[n,1]).*(obj.system.C*xd); % Kx+Cv

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
            s = s.*x(expind(j),:).^ind(expind(j));
        end
        s = repmat(s, [n, 1]);
        y2 = y2+coeff.*s;
    end
    y2 = repmat(p(2,:),[n,1]).*y2;
    y = [y1x; -obj.system.M\(y1d+y2)];
end
end