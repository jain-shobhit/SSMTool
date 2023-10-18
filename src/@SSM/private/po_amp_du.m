function [data, J] = po_amp_du(prob, data, u)
% PO_AMP: This function returns the amplitude of periodic orbits, which is
% defined as the L2 norm of the trajectory (for given degrees-of-freedom).
% Both time-independent (TI) and time-varying (TV) reductions are
% supported.
%
% [DATA, J] = PO_AMP_DU(PROB, DATA, U)
%
% DATA: r,nonauto,coordinates,

% pre-processing
hatr = data.hatr; % hatr.val, hatr.pos (point to pos in c and d), namely <c(pos)-d(pos),r>=hatr.val
cind = full(data.cind);
dind = full(data.dind);
m    = size(cind,2);
coeffs  = data.coeffs;
nhatr   = numel(hatr.val);
optdof = data.optdof;
ndofs   = numel(optdof);
whatr   = zeros(ndofs,nhatr);
isnonauto = data.isnonauto;
isl2norm  = data.isl2norm;
if ~isl2norm
    omega = u(2*m+1); t = u(end);
end

% account for the contribution of nonautonomous part to whatr
if isnonauto
    omega = u(2*m+1); epf = u(2*m+2);
    % solve linear equation to yield x0
    x0 = (1i*data.kappa*omega*data.B-data.A)\data.Flead;
%             x0 = lsqminnorm((data.A-1i*omega*data.B),data.Flead);
    if data.isbaseForce; x0 = x0*omega^2; end
    idplus = find(hatr.val==1);
    whatr(:,idplus) = whatr(:,idplus)+epf*x0(optdof);
    idminus = find(hatr.val==-1);
    whatr(:,idminus) = whatr(:,idminus)+epf*conj(x0(optdof));
end

% Jacobian w.r.t. q=(\rho,\theta) or q=(q^Re,q^Im)
Qbar = data.Qbar;
Jq = zeros(1,2*m);
switch data.coordinates
    case 'polar'
        % u=(\rho,\theta), and (\epsilon,\Omega) if nonauto=true
        rho = u(1:2:2*m-1);
        th  = u(2:2:2*m);
        dAdrho = zeros(m,1);
        dAdth  = zeros(m,1);
        for i=1:nhatr
            idx = hatr.pos{i};
            ci = cind(idx,:); di = dind(idx,:);
            ang = (ci-di)*th;
            rhopower = rho.^((ci+di)');
            pdrho = prod(rhopower,1);
            whatri = coeffs(:,idx).*(pdrho.*exp(1i*ang'));
            whatr(:,i) = whatr(:,i)+sum(whatri,2);
            if isl2norm
                temp = whatr(:,i)'*(Qbar*whatri);
            else
                temp = exp(1i*hatr.val(i)*omega*t)*whatri; 
            end
            dAdrhoi = (temp.*(ci+di)')./rho;
            dAdrho  = dAdrho+sum(dAdrhoi,2);
            dAdthi  = temp.*(ci-di)';
            dAdth   = dAdth+sum(dAdthi,2);
        end  
        dAdth = 1i*dAdth;
        Jq(1:2:end-1) = real(dAdrho);
        Jq(2:2:end)   = real(dAdth);

    case 'cartesian'
        % u=(zRe,zIm), and (\epsilon,\Omega) if nonauto=true
        zRe = u(1:2:2*m-1);
        zIm = u(2:2:2*m);
        zcomp = zRe+1i*zIm; % qs in formulation
        zconj = conj(zcomp);
        zcomp(abs(zcomp)==0) = eps; % handle divided by zero
        zconj(abs(zconj)==0) = eps;
        dAdqRe = zeros(m,1);
        dAdqIm = zeros(m,1);
        for i=1:nhatr
            idx = hatr.pos{i};
            ci = cind(idx,:)'; di = dind(idx,:)';
            zk = prod(zcomp.^ci.*zconj.^di,1);
            whatri = coeffs(:,idx).*zk;
            whatr(:,i) = whatr(:,i)+sum(whatri,2);
            if isl2norm
                temp = whatr(:,i)'*(Qbar*whatri);
            else
                temp = exp(1i*hatr.val(i)*omega*t)*whatri;
            end
            dAdqi      = (temp.*ci)./zcomp;
            dAdqibar   = (temp.*di)./zconj;
            dAdqRe     = dAdqRe+sum(dAdqi+dAdqibar,2);
            dAdqIm     = dAdqIm+sum(dAdqi-dAdqibar,2);            
        end
        dAdqIm = 1i*dAdqIm;
        Jq(1:2:end-1) = real(dAdqRe);
        Jq(2:2:end)   = real(dAdqIm);

    otherwise
        error('Optional coordinates: polar and cartesian');
end

% Jacobian w.r.t omega, epsilon and t
dAdom = []; dAdepf = [];
if ~isl2norm
    hatrval = full(hatr.val');
    dAdt  = sum(1i*omega*whatr.*hatrval.*exp(1i*hatrval*omega*t));
    dAdom = sum(1i*t*whatr.*hatrval.*exp(1i*hatrval*omega*t));
    dAdt  = real(dAdt);
    dAdom = real(dAdom);
end
if isnonauto
    dom = (data.A-1i*data.kappa*omega*data.B)\(data.B*x0);
    if isl2norm
        dAdom  = 2*epf*real(1i*whatr(:,idplus)'*(Qbar*dom(optdof)));
        if data.isbaseForce
            dAdom = dAdom+4*epf*real(whatr(:,idplus)'*(Qbar*x0(optdof)))/omega;
        end
        dAdepf = 2*real(whatr(:,idplus)'*(Qbar*x0(optdof)));
    else
        dAdom  = dAdom+2*epf*real(1i*dom(optdof)*exp(1i*omega*t));
        if data.isbaseForce
            dAdom = dAdom+4*epf*real(x0(optdof)*exp(1i*omega*t))/omega;
        end
        dAdepf = 2*real(x0(optdof)*exp(1i*omega*t));
    end
end

% calculate response amplitude
if isl2norm
    A = sqrt(real(sum(sum(conj(whatr).*(data.Q*whatr)))));
    J = [Jq dAdom dAdepf]/A;
else
    J = [Jq dAdom dAdepf dAdt];
end
J = real(J);

% tic
% [data, Jd] = coco_ezDFDX('f(o,d,x)', prob, data, @po_amp, u);
% toc
% max(max(abs(J-Jd(:)')))

end