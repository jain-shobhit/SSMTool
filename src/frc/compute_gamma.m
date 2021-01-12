function gamma = compute_gamma(R_0)
order = length(R_0);

n_gamma = floor((order-1)/2);
gamma = zeros(n_gamma,1);
for j = 1:n_gamma
    if ~isempty(R_0{2*j+1}.ind)
        [~, loc] = ismember([j+1,j],R_0{2*j+1}.ind,'rows');
        gamma(j) = R_0{2*j+1}.coeffs(1,loc);
    end
end
fprintf('gamma = \n')
disp(gamma)