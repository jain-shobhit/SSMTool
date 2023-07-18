function fext = compute_fext(obj,t,x,xd)
% COMPUTE_FEXT We compute the external force at a particular time t 
% in a second-order mechanical system. 
% COMPUTE_FEXT We compute the external force at a particular time t
% in a second-order mechanical system.

% ﻿The quasiperiodic second order system force is given as a Taylor expansion on coordinates of
%  the physical coordinates:
%   $ \mathbf {f}(\mathbf{x},\text{\boldmath$\phi$}) = \begin{bmatrix} f^1(\mathbf{x},\text{\boldmath$\phi$})
%     \\ \vdots \\  f^{n}(\mathbf{x},\text{\boldmath$\phi$})\end{bmatrix}, \ f^i(\mathbf{x},\text{\boldmath$\phi$})
%     = \sum_{\mathbf{n}\in \mathbb{N}^{n}}
%     f^i_{\mathbf{n}}(\text{\boldmath$\phi$}) \mathbf{x}^\mathbf{n} $
% ﻿The force coefficients are given as a fourier expansion in terms of the phase variable $\text{\boldmath$\phi$ = \Omegat}$
%   $ f^b_{\mathbf{k}}(\text{\boldmath$\phi$}) = \sum_{\text{\boldmath$\eta$}
%  \in \mathbb{Z}^k} f^b_{\mathbf{k},\text{\boldmath$\eta$} } e^{i\langle \text{\boldmath$\eta$},
%   \text{\boldmath$\phi$}\rangle} $


assert(obj.order == 2, ' fext can only be computed for second-order systems')

if isempty(obj.fext)
    fext = sparse(obj.n,1);
else
    assert(~isempty(obj.Omega), ' fext cannot be evaluated as the Omega property of the DS class is empty')
    % This function assumes periodic forcing
    nt = size(x,2);

    fext = zeros(obj.n,nt);

    nKappa = obj.nKappa;

    for i = 1:nKappa
        % Highest order this Force contributes at
        order = numel(obj.fext.data(i).f_n_k)-1;

        % zeroth order
        f0 = obj.fext.data(i).f_n_k(1);

        if ~isempty(f0) && ~isempty(f0.coeffs)            

            fext = fext + real( f0.coeffs * exp(1i  * obj.fext.data(i).kappa * obj.Omega .* t));
        end

        for j = 1:order % array starting at 0
            %Contribution to order j Forcing with harmonic kappa_i
            fij = obj.fext.data(i).f_n_k(j+1);
            if ~isempty(fij) && ~isempty(fij.coeffs) 
                if size(fij.ind,2) == obj.n
                    fext = fext + real(expand_multiindex(fij,x) .* exp(1i  * obj.fext.data(i).kappa  * obj.Omega .* t));
                else
                    fext = fext + real(expand_multiindex(fij,[x;xd]) .* exp(1i  * obj.fext.data(i).kappa  * obj.Omega .* t));
                end

            end
        end
    end
    assert(size(obj.Omega,2) == size(obj.fext.epsilon,2),'epsilon and Omega must have same number of columns'); 
    fext = fext * diag(obj.fext.epsilon);
    if obj.Options.BaseExcitation
        fext = fext*diag((obj.Omega).^2);
    end
end

