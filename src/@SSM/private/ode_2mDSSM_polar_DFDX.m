function J = ode_2mDSSM_polar_DFDX(z, p, data)
% ODE_2MDSSM_POLAR_DFDX This function presents vectorized implementation of
% the Jacobian of the vector field of the reduced dynamics on
% 2m-dimensional SSMs with respect to state z. Here z is a 
% 2m-dimensinoal state vector and p is parameter vector for excitation
% frequency and amplitude. All other info such as eigenvalues of master
% subspace, coefficients of nonlinear terms is included in the structure
% data. The state vector here is in the form of polar coordinate.
%
% See also: ODE_2MDSSM_CARTESIAN_DFDX
% tic
% extract data fields
beta   = data.beta;
kappa  = data.kappa;
lamdRe = data.lamdRe;
lamdIm = data.lamdIm;
mFreqs = data.mFreqs;
iNonauto = data.iNonauto;
rNonauto = data.rNonauto;
rNonauto = full(rNonauto);

% rename state and parameter
rho = z(1:2:end-1,:);
th  = z(2:2:end,:);
om   = p(1,:);
epsf = p(2,:);

m  = numel(mFreqs);
nt = numel(om);
J  = zeros(2*m,2*m,nt); 
% autonomous part
em = eye(m);
for i=1:m
    % linear part
    J(2*i-1,2*i-1,:) = lamdRe(i);
    % nonlinear part
    kappai = kappa{i};
    kappai = full(kappai);
    betai  = beta{i};
    nka = size(kappai,1);
    ei = em(i,:);
    for k=1:nka
        ka = kappai(k,:);
        be = betai(k);
        l = ka(1:2:end-1);
        j = ka(2:2:end);
        ang = (l-j-ei)*th;
        rhopower = rho.^((l+j)');
        pdrho = prod(rhopower,1);
        drho1 = (l+j)'./rho;
        % df_rho/drho
        J(2*i-1,1:2:end-1,:) = J(2*i-1,1:2:end-1,:)+...
            reshape(drho1.*(pdrho.*(real(be)*cos(ang)-imag(be)*sin(ang))),[1,m,nt]);
        % df_rho/dtheta
        J(2*i-1,2:2:end,:) = J(2*i-1,2:2:end,:)+...
            reshape((l-j-ei)'.*(pdrho.*(-real(be)*sin(ang)-imag(be)*cos(ang))),[1,m,nt]);
        % df_theta/drho
        drho2 = (l+j-ei)'./rho;
        J(2*i,1:2:end-1,:) = J(2*i,1:2:end-1,:)+...
            reshape(drho2.*(pdrho./rho(i,:).*(real(be)*sin(ang)+imag(be)*cos(ang))),[1,m,nt]);
        % df_theta/dtheta
        J(2*i,2:2:end,:) = J(2*i,2:2:end,:)+...
            reshape((l-j-ei)'.*(pdrho./rho(i,:).*(real(be)*cos(ang)-imag(be)*sin(ang))),[1,m,nt]);
    end
end

% nonautonomous leading part
for i=1:numel(iNonauto)
    id = iNonauto(i);
    r  = rNonauto(i);
    r  = epsf*r; 
    if data.isbaseForce; r = r.*om.^2; end
    rRe = real(r);
    rIm = imag(r);
    J(2*id-1,2*id,:) = J(2*id-1,2*id,:)+reshape(-rRe.*sin(th(id,:))+rIm.*cos(th(id,:)),[1,1,nt]);
    J(2*id,2*id-1,:) = J(2*id,2*id-1,:)+reshape(-(-rRe.*sin(th(id,:))+rIm.*cos(th(id,:)))./rho(id,:).^2,[1,1,nt]);
    J(2*id,2*id,:)   = J(2*id,2*id,:)+reshape(-(rRe.*cos(th(id,:))+rIm.*sin(th(id,:)))./rho(id,:),[1,1,nt]);
end

% func = @(z,p) ode_2mDSSM_polar(z,p,data);
% JJ = coco_ezDFDX('f(x,p)v',func,z,p);
% ers = max(abs(J(:)-JJ(:)))

end