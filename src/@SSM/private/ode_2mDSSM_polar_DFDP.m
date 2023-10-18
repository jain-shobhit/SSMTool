function J = ode_2mDSSM_polar_DFDP(z, p, data)
% ODE_2MDSSM_POLAR_DFDX This function presents vectorized implementation of
% the Jacobian of the vector field of the reduced dynamics on
% 2m-dimensional SSMs with respect to system parameters p. Here z is a 
% 2m-dimensinoal state vector and p is parameter vector for excitation
% frequency and amplitude. All other info such as eigenvalues of master
% subspace, coefficients of nonlinear terms is included in the structure
% data. The state vector here is in the form of polar coordinate.
%
% See also: ODE_2MDSSM_CARTESIAN_DFDX

% extract data fields
mFreqs = data.mFreqs;
iNonauto = data.iNonauto;
rNonauto = data.rNonauto;

% rename state and parameter
rho = z(1:2:end-1,:);
th  = z(2:2:end,:);
om  = p(1,:);
nt  = numel(om); 
epsf = p(2,:);

J = zeros(size(z,1),2,nt);

% w.r.t. Omega
J(2:2:end,1,:) = repmat(-mFreqs(:),[1,1,nt]);
% w.r.t. epsilon
m = numel(mFreqs);
yrho = zeros(m,nt);
yth  = zeros(m,nt);
for i=1:numel(iNonauto)
    id = iNonauto(i);
    r  = rNonauto(i);
    rRe = real(r);
    rIm = imag(r);
    yrho(id,:) = yrho(id,:)+rRe.*cos(th(id,:))+rIm.*sin(th(id,:));
    yth(id,:)  = yth(id,:)-rRe.*sin(th(id,:))./rho(id,:)+rIm.*cos(th(id,:))./rho(id,:);
end
J(1:2:end-1,2,:) = yrho;
J(2:2:end,2,:)   = yth;

if data.isbaseForce
    J(1:2:end-1,1,:) = J(1:2:end-1,1,:)+reshape(2*epsf.*om.*yrho,[m,1,nt]);
    J(2:2:end,1,:)   = J(2:2:end,1,:)+reshape(2*epsf.*om.*yth,[m,1,nt]);
    J(1:2:end-1,2,:) = J(1:2:end-1,2,:).*reshape(repmat(om.^2,[m,1]),[m,1,nt]);
    J(2:2:end,2,:)   = J(2:2:end,2,:).*reshape(repmat(om.^2,[m,1]),[m,1,nt]);
end


% func = @(z,p) ode_2mDSSM_polar(z,p,data);
% JJ = coco_ezDFDP('f(x,p)v',func,z,p);
% ers = max(abs(J(:)-JJ(:)))
end