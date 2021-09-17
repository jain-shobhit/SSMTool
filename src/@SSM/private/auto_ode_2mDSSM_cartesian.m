function y = auto_ode_2mDSSM_cartesian(z, data)
% ODE_2MDSSM_CARTESIAN This function presents vectorized implementation of vector
% field of reduced dynamics on 2m-dimensional SSMs. Here z is a
% 2m-dimensional state vector and p is parameter vector for excitation
% frequency and amplitude. All other info such as eigenvalues of master
% spectral subspace, coefficients of nonlinear terms is included in the
% structure data. The state vector here is in the form of Cartesian
% coordinates.
%
% See also: ODE_2MDSSM_POLAR

assert(~isempty(data), 'Structure data in ode_2mDSSM_cartesian is empty');
% extract data fields
beta   = data.beta;
kappa  = data.kappa;
lamdRe = data.lamdRe;
lamdIm = data.lamdIm;

% rename state and parameter
zRe  = z(1:2:end-1,:);
zIm  = z(2:2:end,:);
zcomp = zRe+1i*zIm; % qs in formulation
zconj = conj(zcomp);

% autonomous linear part
yRe = lamdRe.*zRe-lamdIm.*zIm;
yIm = lamdRe.*zIm+lamdIm.*zRe;

% autonomous nonlinear part
m  = numel(lamdRe);
for i=1:m
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
        yRe(i,:) = yRe(i,:)+real(zk);
        yIm(i,:) = yIm(i,:)+imag(zk);
    end
end

nt = size(z,2);
y  = zeros(2*m,nt);
y(1:2:end-1,:) = yRe;
y(2:2:end,:)   = yIm;

end