function J = ode_2mDSSM_cartesian_DFDP(z, p, data)
% ODE_2MDSSM_POLAR_DFDX This function presents vectorized implementation of
% the Jacobian of the vector field of the reduced dynamics on
% 2m-dimensional SSMs with respect to system parameters p. Here z is a 
% 2m-dimensinoal state vector and p is parameter vector for excitation
% frequency and amplitude. All other info such as eigenvalues of master
% subspace, coefficients of nonlinear terms is included in the structure
% data. The state vector here is in the form of polar coordinate.
%
% See also: ODE_2MDSSM_POLAR_DFDX

% extract data fields
mFreqs = data.mFreqs;
mFreqs = mFreqs(:);
iNonauto = data.iNonauto;
rNonauto = data.rNonauto;

% rename state and parameter
zRe = z(1:2:end-1,:);
zIm = z(2:2:end,:);
om  = p(1,:);
epsf= p(2,:);
nt  = numel(om); 

J = zeros(size(z,1),2,nt);

% w.r.t. Omega
J(1:2:end-1,1,:) = mFreqs.*zIm;
J(2:2:end,1,:)   = -mFreqs.*zRe;
% w.r.t. epsilon
% nonautonomous leading part
m = numel(mFreqs);
yRe = zeros(m,nt);
yIm = zeros(m,nt);
for i=1:numel(iNonauto)
    id = iNonauto(i);
    r  = rNonauto(i);
    rRe = real(r);
    rIm = imag(r);
    yRe(id,:) = yRe(id,:)+rRe;
    yIm(id,:) = yIm(id,:)+rIm;
end
J(1:2:end-1,2,:) = yRe;
J(2:2:end,2,:)   = yIm;

if data.isbaseForce
    J(1:2:end-1,1,:) = J(1:2:end-1,1,:)+reshape(2*epsf.*om.*yRe,[m,1,nt]);
    J(2:2:end,1,:)   = J(2:2:end,1,:)+reshape(2*epsf.*om.*yIm,[m,1,nt]);
    J(1:2:end-1,2,:) = J(1:2:end-1,2,:).*reshape(repmat(om.^2,[m,1]),[m,1,nt]);
    J(2:2:end,2,:)   = J(2:2:end,2,:).*reshape(repmat(om.^2,[m,1]),[m,1,nt]);
end

% func = @(z,p) ode_2mDSSM_cartesian(z,p,data);
% JJ = coco_ezDFDP('f(x,p)v',func,z,p);
% ers = max(abs(J(:)-JJ(:)))
end