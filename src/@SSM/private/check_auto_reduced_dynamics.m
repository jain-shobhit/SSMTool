function [beta,kappa] = check_auto_reduced_dynamics(R_0,order,mFreqs)
% check reduced dynamics (factor out exp(i*Omega*t) part
m  = numel(mFreqs);
Em = eye(m);
beta  = cell(m,1); % coefficients - each cell corresponds to one mode
kappa = cell(m,1); % exponants
for k = 2:order
    R = R_0{k};
    coeffs = R.coeffs;
    ind = R.ind;
    if ~isempty(coeffs)
        for i=1:m
            betai = coeffs(2*i-1,:);
            [~,ki,betai] = find(betai);
            kappai = ind(ki,:);
            % check resonant condition  
            l = kappai(:,1:2:end-1);
            j = kappai(:,2:2:end);
            nk = numel(ki);
            rm = repmat(mFreqs(:)',[nk,1]);
            flagi = dot(l-j-repmat(Em(i,:),[nk,1]),rm,2);
            assert(all(flagi==0), 'Reduced dynamics is not consisent with desired IRs');
            % assemble terms
            beta{i}  = [beta{i} betai];
            kappa{i} = [kappa{i}; kappai];
        end
    end
end
end
