function J = ode_2mDSSM_cartesian_DFDX(z, p, data)
% ODE_2MDSSM_POLAR_DFDX This function presents vectorized implementation of
% the Jacobian of the vector field of the reduced dynamics on
% 2m-dimensional SSMs with respect to state z. Here z is a 
% 2m-dimensinoal state vector and p is parameter vector for excitation
% frequency and amplitude. All other info such as eigenvalues of master
% subspace, coefficients of nonlinear terms is included in the structure
% data. The state vector here is in the form of polar coordinate.
%
% See also: ODE_2MDSSM_POLAR_DFDX

% extract data fields
beta   = data.beta;
kappa  = data.kappa;
lamdRe = data.lamdRe;
lamdIm = data.lamdIm;
mFreqs = data.mFreqs;
iNonauto = data.iNonauto;
rNonauto = data.rNonauto;

% rename state and parameter
zRe  = z(1:2:end-1,:);
zIm  = z(2:2:end,:);
om   = p(1,:);
epsf = p(2,:);
zcomp = zRe+1i*zIm; % qs in formulation
zconj = conj(zcomp);
zcomp(zcomp==0) = eps; % to handle 0/0 later
zconj(zconj==0) = eps;

m  = numel(mFreqs);
nt = numel(om);
J  = zeros(2*m,2*m,nt); 
% autonomous part

% autonomous linear part
mom = mFreqs(:).*om;
yRe = lamdRe.*zRe-lamdIm.*zIm+zIm.*mom;
yIm = lamdRe.*zIm+lamdIm.*zRe-zRe.*mom;

% autonomous nonlinear part
m  = numel(mFreqs);
for i=1:m
    % linear part
    J(2*i-1,2*i-1,:) = lamdRe(i);
    J(2*i-1,2*i,:)   = mFreqs(i)*om-lamdIm(i);
    J(2*i,2*i-1,:)   = lamdIm(i)-mFreqs(i)*om;
    J(2*i,2*i,:)     = lamdRe(i);
    % nonlinear part
    kappai = kappa{i};
    kappai = full(kappai);
    betai  = beta{i};
    nka = size(kappai,1);
    nbe = numel(betai);
    assert(nka==nbe, 'Size of kappa%d and beta%d does not match',i,i);
    for k=1:nka
        ka = kappai(k,:);
        be = betai(k);
        l = ka(1:2:end-1)';
        j = ka(2:2:end)';
        zk = be*prod(zcomp.^l.*zconj.^j,1);
        dz1 = l./zcomp+j./zconj; 
        dz2 = l./zcomp-j./zconj;
        % df_qR/dqR
        J(2*i-1,1:2:end-1,:) = J(2*i-1,1:2:end-1,:)+reshape(real(dz1.*zk),[1,m,nt]);
        % df_qR/dqI
        J(2*i-1,2:2:end,:)   = J(2*i-1,2:2:end,:)+reshape(real(1j*dz2.*zk),[1,m,nt]);
        % df_qI/dqR
        J(2*i,1:2:end-1,:)   = J(2*i,1:2:end-1,:)+reshape(imag(dz1.*zk),[1,m,nt]);
        % df_qI/dqI
        J(2*i,2:2:end,:)     = J(2*i,2:2:end,:)+reshape(imag(1j*dz2.*zk),[1,m,nt]);      
    end
end

% nonautonomous leading part - no contributions to Jacobian
% 
% func = @(z,p) ode_2mDSSM_cartesian(z,p,data);
% JJ = coco_ezDFDX('f(x,p)v',func,z,p);
% ers = max(abs(J(:)-JJ(:)))
end