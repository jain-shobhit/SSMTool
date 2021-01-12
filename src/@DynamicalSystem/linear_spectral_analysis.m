function [V, D, W] = linear_spectral_analysis(obj)

if obj.N < obj.Options.Nmax
    if ~issparse(obj.A)
        [V, LAMBDA, W] = eig(obj.A,obj.B);
    else
        [V, LAMBDA, W] = eig(full(obj.A),full(obj.B));
    end
else
    
    E_max = obj.Options.Emax;
    n = obj.n;
    
    if obj.order == 2
        
    disp(['Due to high-dimensionality, we compute only the first ' num2str(E_max) ' eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients'] )
    % Computing undamped eigenvalues and eigenvectors
    
    disp ('Assuming a proportional damping hypthesis with symmetric matrices') 
    [U, ~, NOT_CONVERGED] = eigs(sparse(obj.K),sparse(obj.M),E_max,'smallestabs');
    
    % Assuming proportional damping for estimating damped eigenvectors
    LAMBDA = zeros(2*E_max,2*E_max);
    V = zeros(2*n,2*E_max);
    
    % Rearranging eigenvalues
    for j = 1:E_max
        mu_j = U(:,j).'* obj.M * U(:,j);
        omega2_j = (U(:,j).'* obj.K * U(:,j))/mu_j;
        beta_j = (U(:,j).'* obj.C * U(:,j))/mu_j;
        
        fprintf('modal damping ratio for %d mode is %d\n', j, beta_j/(2*sqrt(omega2_j)));
        
        lambda1 = (-beta_j + sqrt(beta_j^2 - 4 * omega2_j) ) / 2;
        lambda2 = (-beta_j - sqrt(beta_j^2 - 4 * omega2_j) ) / 2;
        
                
        LAMBDA(2*j-1,2*j-1) = lambda1;
        LAMBDA(2*j,2*j) = lambda2;
        
        V(:,2*j-1) = [U(:,j); lambda1*U(:,j)];
        V(:,2*j) = [U(:,j); lambda2*U(:,j)];
    end
    W = conj(V);
    
    if NOT_CONVERGED
        error('The eigenvalue computation did not converge, please adjust the number of eigenvalues to be computed')
    end
    
    if ~issymmetric(obj.K) || ~issymmetric(obj.M)
        disp('the left eigenvectors may be incorrect in case of asymmetry of matrices')
%         [W, ~] = eigs(A',B',E_max,'smallestabs');
        % we assume here that both eigenvalue problems return eigenvalues in
        % the same order.
    end
    
    else
        % right eigenvectors
        [V, Dv] = eigs(obj.A,obj.B,E_max,'smallestabs');
        [Lambda_sorted,I] = sort(diag(Dv),'descend','ComparisonMethod','real');
        LAMBDA = diag(Lambda_sorted);
        V = V(:,I);
        
        % left eigenvectors
        if issymmetric(obj,A) && issymmetric(obj.B)
            W = conj(V);
        else
            [W, Dw] = eigs(obj.A',obj.B',E_max,'smallestabs');
            [~,I] = sort(diag(Dw),'descend','ComparisonMethod','real');
            W = W(:,I);
            W = conj(W);
        end
    end
    
end

[V, D, W] = sort_modes(V, LAMBDA, W);
[V,W] =  normalize_modes(V,W,obj.B);

obj.spectrum.V = V;
obj.spectrum.W = W;
obj.spectrum.Lambda = D;

fprintf('\n The first %d eigenvalues are given as \n',length(D))
disp(D)

end

function [V, Lambda, W] = sort_modes(V, D, W)
%SORT_MODES: This function sorts the eigenvectors (V) in descending order of
%the real parts of the corresponding eigenvalues (D). The resulting
%eigenvectors are also normalized to unit magnitude. 

% obtain the eigenvalues as a vector instead of a diagonal matrix
Lambda = diag(D); 
if ~iscolumn(Lambda)
    Lambda = transpose(Lambda);
end
% sort eigenvalues in the descending order of real parts, incase of tie by
% ascending order of magnitude of imaginary parts
[Lambda_sorted,I] = sortrows([real(Lambda), abs(imag(Lambda)) sign(imag(Lambda))],[1 2],{'descend' 'ascend'});
D = Lambda_sorted(:,1) + 1i * Lambda_sorted(:,2).*Lambda_sorted(:,3);
D = diag(D);
% arrange eigenvectors accordingly
V = V(:,I);
W = W(:,I);

% ensure positive imaginary part first in every complex pair
Lambda = diag(D);
skip = false;
for j = 1:length(Lambda)
    if skip 
        skip = false;
        continue;
    end    
    if ~isreal(Lambda(j))
        % extract complex eigenpair
        V0 = V(:,j:j+1);
        W0 = W(:,j:j+1);
        Lambda0 = Lambda(j:j+1);
        % sort eigenvalues in the descending order of imaginary parts
        [~,I] = sort(imag(Lambda0),'descend','ComparisonMethod','real');        
        Lambda(j) = Lambda0(I(1));
        V(:,j) = V0(:,I(1));
        W(:,j) = W0(:,I(1));
        % ensure complex conjugate eigenvalues and eigenvectors 
        Lambda(j+1) = conj(Lambda(j));                
        V(:,j+1) = conj(V(:,j));
        W(:,j+1) = conj(W(:,j));
        % move to the next pair of eigenvalues
        skip = true; 
    end
end
% D = diag(Lambda);
end

function [V,W] = normalize_modes(V,W,B)
%NORMALIZE_MODES: This function normalizes the right and left eigenvectors  
% W, V withrespect to the matrix B. 

% V = V*diag(1./vecnorm(V));
mu = diag(W'*B*V);

% V = V*diag(1./sqrt(mu));
% W = W*diag(1./(sqrt(mu)'));

W = W*diag(1./(mu'));
end
