function [data, y] = po_amp(prob, data, u)
% PO_AMP: This function returns the amplitude of periodic orbits, which is
% defined as the L2 norm of the trajectory (for given degrees-of-freedom).
% Both time-independent (TI) and time-varying (TV) reductions are
% supported.
%
% [DATA, Y] = PO_AMP(PROB, DATA, U)
%
% DATA: r,nonauto,coordinates,

% pre-processing
hatr = data.hatr; % hatr.val, hatr.pos (point to pos in c and d), namely <c(pos)-d(pos),r>=hatr.val
cind = data.cind;
dind = data.dind;
m    = size(cind,2);
coeffs  = data.coeffs;
nhatr   = numel(hatr.val);
optdof = data.optdof;
ndofs   = numel(optdof);
whatr   = zeros(ndofs,nhatr);
isnonauto = data.isnonauto;

% autonomous part
switch data.coordinates
    case 'polar'
        % u=(\rho,\theta), and (\epsilon,\Omega) if nonauto=true
        rho = u(1:2:2*m-1);
        th  = u(2:2:2*m);
        for i=1:nhatr
            idx = hatr.pos{i};
            ci = cind(idx,:); di = dind(idx,:);
            ang = (ci-di)*th;
            rhopower = rho.^((ci+di)');
            pdrho = prod(rhopower,1);
            whatri = coeffs(:,idx).*(pdrho.*exp(1i*ang'));
            whatr(:,i) = sum(whatri,2);
        end       
        
    case 'cartesian'
        % u=(zRe,zIm), and (\epsilon,\Omega) if nonauto=true
        zRe = u(1:2:2*m-1);
        zIm = u(2:2:2*m);
        zcomp = zRe+1i*zIm; % qs in formulation
        zconj = conj(zcomp);
        for i=1:nhatr
            idx = hatr.pos{i};
            ci = cind(idx,:)'; di = dind(idx,:)';
            zk = prod(zcomp.^ci.*zconj.^di,1);
            whatri = coeffs(:,idx).*zk;
            whatr(:,i) = sum(whatri,2);
        end
        
    otherwise
        error('Optional coordinates: polar and cartesian');
end

% nonautonomous part
if isnonauto
    omega = u(2*m+1); epf = u(2*m+2);
    % solve linear equation to yield x0
    x0 = (1i*data.kappa*omega*data.B-data.A)\data.Flead;
%             x0 = lsqminnorm((1i*data.kappa*omega*data.B-data.A),data.Flead,TOL);
    if data.isbaseForce; x0 = x0.*omega.^2; end
    idplus = find(hatr.val==1);
    whatr(:,idplus) = whatr(:,idplus)+epf*x0(optdof);
    idminus = find(hatr.val==-1);
    whatr(:,idminus) = whatr(:,idminus)+epf*conj(x0(optdof));
end
% calculate response amplitude
if data.isl2norm
    y = sqrt(real(sum(sum(conj(whatr).*(data.Q*whatr)))));
else
    omega = u(2*m+1);
    t = u(end);
    hatrval = full(hatr.val');
    y = exp(1i*hatrval*omega*t);
    y = sum(whatr.*y);
    y = real(y);
end

end