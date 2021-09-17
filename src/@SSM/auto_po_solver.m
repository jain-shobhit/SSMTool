function varargout = auto_po_solver(obj,R_0,oid,t0,p0,coordinates)
% AUTO_PO_SOLVER This function solves the periodic orbit of autonomous
% dynamics on SSM with given initial guess of periodic orbit {t0,p0}. Here
% t0 is a time sequence, and p0 are the corresponding sequence of states in
% either polar or cartesian representation (string in coordinates). For the
% rest of arguments, lambda is the spectrum of the spectral subspace used
% for yielding the reduced dynamics R_0

lambda = obj.System.spectrum.Lambda;
lambdaRe = real(lambda);
lambdaIm = imag(lambda);

% check reduced dynamics to see its consistent with reduced dynamics
m     = numel(lambda)/2;
order = numel(R_0);
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
            % assemble terms
            beta{i}  = [beta{i} betai];
            kappa{i} = [kappa{i}; kappai];
        end
    end
end

fdata = struct();
fdata.beta  = beta;
fdata.kappa = kappa;
fdata.lamdRe = lambdaRe(1:2:end-1);
fdata.lamdIm = lambdaIm(1:2:end-1);

ispolar = strcmp(coordinates, 'polar');
fdata.ispolar = ispolar;
if ispolar
    error('autonomous vector field in polar representation is not implemented yet');
else
    odefun = @(z,p) auto_ode_2mDSSM_cartesian(z,fdata);
end
funcs  = {odefun};

prob = coco_prob();
prob = coco_set(prob,'corr','TOL',1e-5);
prob = ode_isol2po(prob, '', funcs{:}, t0, p0);
runid = coco_get_id(oid, 'auto-po');
bd = coco(prob, runid, [], 0);
varargout{1} = bd;
end


