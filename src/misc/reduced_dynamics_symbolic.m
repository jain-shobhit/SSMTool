function y = reduced_dynamics_symbolic(lamdMaster,R0,type)
%REDUCED_DYNAMICS_SYMBOLIC This function returns vector field in symbolic
%expression, which can be used for latex documentation. In addition, it
%write the vector field to a function file in hard disk.
%
% y = reduced_dynamics_symbolic(lamdMaster,R0,type)
%
% lamdMaster - spectrum of master spectral subspace
% R0 - Reduced dynamics from compute_whisker
% type - a string, which should be either 'auto' or 'nonauto'

lamdRe = real(lamdMaster);
lamdIm = imag(lamdMaster);
lamdRe = lamdRe(1:2:end-1);
lamdIm = lamdIm(1:2:end-1);
order  = numel(R0);
m      = numel(lamdRe);
beta   = cell(m,1); % coefficients - each cell corresponds to one mode
kappa  = cell(m,1); % exponants

% create symbolic array
rho = []; theta = [];
for k=1:m
    rho = [rho; sym(['rho_',num2str(k)])];
    theta = [theta; sym(['theta_',num2str(k)])];
end

switch type
    case 'auto'
        % assemble and simplify nonlinear coefficients
        for k = 2:order
            R = R0{k};
            coeffs = R.coeffs;
            ind = R.ind;
            if ~isempty(coeffs)
                for i=1:m
                    betai = coeffs(2*i-1,:);
                    % remove small coefficients
                    Rebetai = real(betai);
                    Imbetai = imag(betai);
                    bounds  = max(1e-6*norm(betai),1e-8);
                    Rebetai(abs(Rebetai)<bounds) = 0;
                    Imbetai(abs(Imbetai)<bounds) = 0;
                    betai = Rebetai+1j*Imbetai;
                    [~,ki,betai] = find(betai);
                    kappai = ind(ki,:);
                    % assemble terms
                    beta{i}  = [beta{i} betai];
                    kappa{i} = [kappa{i}; kappai];
                end
            end
        end
        % formualtion of vector field
        % remove small parts of eigenvalues
        lamdRe(abs(lamdRe)<1e-6*norm(lamdMaster)) = 0;
        lamdIm(abs(lamdIm)<1e-6*norm(lamdMaster)) = 0;
        yrho = lamdRe.*rho;
        yth  = sym(lamdIm);
        em = eye(m);
        for i=1:m
            kappai = kappa{i};
            kappai = full(kappai);
            betai  = beta{i};
            nka = size(kappai,1);
            nbe = numel(betai);
            assert(nka==nbe, 'Size of kappa%d and beta%d does not match',i,i);
            ei = em(i,:);
            for k=1:nka
                ka = kappai(k,:);
                be = betai(k);
                l = ka(1:2:end-1);
                j = ka(2:2:end);
                ang = (l-j-ei)*theta;
                rhopower = rho.^((l+j)');
                pdrho = prod(rhopower,1);
                yrho(i) = yrho(i)+pdrho*(real(be)*cos(ang)-imag(be)*sin(ang));
                yth(i)  = yth(i)+pdrho/rho(i)*(real(be)*sin(ang)+imag(be)*cos(ang));
            end
        end
        % write to y
        y = [yrho;yth];
        y(1:2:end-1) = yrho;
        y(2:2:end)   = yth;       
    case 'nonauto'
        error('nonauto is not supported yet');
        
    otherwise
        disp('the second argument should be auto/nonauto in string');
end

end